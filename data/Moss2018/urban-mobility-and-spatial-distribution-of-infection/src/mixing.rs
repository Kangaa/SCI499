use csv;
use ndarray::prelude::*;
use ndarray::Zip;
use std::collections::HashMap;
use std::convert::From;
use std::fs::{File, read_dir};
use std::path::Path;

pub struct OdMat {
    pub m: Array2<f64>,
    pub names: Vec<String>,
    pub popns: Array1<usize>,
}

impl OdMat {
    pub fn read<P>(mm_path: P, popn_path: P) -> OdMat where P: AsRef<Path> {
        let mm_file = File::open(mm_path).unwrap();
        let mut reader = csv::Reader::from_reader(mm_file);
        let mut names: Vec<String> = vec![];
        let mut values: Vec<f64> = vec![];
        for record in reader.records() {
            let row = record.unwrap();
            names.push(row[0].to_string());
            for val in row.iter().skip(1) {
                values.push(val.parse().unwrap());
            }
        }
        let n = names.len();
        let mat = Array::from_shape_vec((n, n), values).unwrap();
        let popns = read_popns(popn_path, &names);
        OdMat {
            m: mat,
            names: names,
            popns: popns,
        }
    }
}

/// Each region can either share the same frac_self value (`Uniform`) or have
/// independent values (`Variable`).
#[derive(Clone, Debug)]
pub enum FracSelf {
    Uniform(f64),
    Variable(Array1<f64>),
}

pub struct MixMat {
    pub m: Array2<f64>,
    pub names: Vec<String>,
    pub popns: Array1<usize>,
    pub frac_self: FracSelf,
    pub frac_cbd:f64,
}

impl MixMat {
    fn scale_diag_var(mat: &mut Array2<f64>, frac_self: &Array1<f64>) {
        Zip::from(mat.genrows_mut())
            .and(frac_self)
            .apply(|mut row, &fs| {
                row *= (1.0 - fs) / row.scalar_sum();
            });
        mat.diag_mut().assign(frac_self);
    }

    fn scale_diag_unif(mat: &mut Array2<f64>, frac_self: f64) {
        let frac_mix = 1.0 - frac_self;
        for mut row in mat.genrows_mut() {
            let scale: f64 = frac_mix / row.scalar_sum();
            row *= scale;
        }
        mat.diag_mut().fill(frac_self);
    }

    fn scale_diag(mat: &mut Array2<f64>, frac_self: &FracSelf) {
        match frac_self {
            &FracSelf::Uniform(val) => Self::scale_diag_unif(mat, val),
            &FracSelf::Variable(ref arr) => Self::scale_diag_var(mat, &arr),
        }
    }

    pub fn test_scale_diag(mat: &Array2<f64>, frac_self: f64) {
        let mut m1 = mat.clone();
        Self::scale_diag_unif(&mut m1, frac_self);

        let mut m2 = mat.clone();
        let mut fs2: Array1<f64> = Array::zeros(m2.rows());
        fs2.fill(frac_self);
        Self::scale_diag_var(&mut m2, &fs2);

        let mut m3 = mat.clone();
        let mut fs3: Array1<f64> = Array::zeros(m3.rows());
        fs3.fill(frac_self);
        fs3[1] *= 0.9;
        Self::scale_diag_var(&mut m3, &fs3);

        for row in m1.genrows() {
            assert!((1.0 - row.scalar_sum()).abs() < 1e-8, "m1 rowsum");
        }
        for row in m2.genrows() {
            assert!((1.0 - row.scalar_sum()).abs() < 1e-8, "m2 rowsum");
        }
        for row in m3.genrows() {
            assert!((1.0 - row.scalar_sum()).abs() < 1e-8, "m3 rowsum");
        }
        assert!(m1 == m2, "mixing matrices differ");
        assert!(m1 != m3, "mixing matrices do not differ");
    }

    /// Handle the CBD differently.
    fn redistribute(mat: &mut Array2<f64>, popns: &Array1<usize>,
                    frac_cbd: f64, cbd_ix: usize) {
        // Determine how much all other regions mix with the CBD.
        let mut mix_w_cbd = mat.column(cbd_ix).to_owned();
        mix_w_cbd[cbd_ix] = 0.0;
        // Weight the mixing by each region's resident population.
        mix_w_cbd *= &popns.map(|&x| x as f64);
        // Normalise these values to conserve the force of infection.
        mix_w_cbd /= mix_w_cbd.scalar_sum();
        for rix in 0..mat.rows() {
            if rix == cbd_ix {
                continue;
            }
            let mut row = mat.row_mut(rix);
            let cbd_mix = row[cbd_ix];
            let add_mix: Vec<f64> = mix_w_cbd.iter()
                .map(|v| cbd_mix * (1.0 - frac_cbd) * v)
                .collect();
            let add_mix: Array1<f64> = From::from(add_mix);
            row += &add_mix;
            row[cbd_ix] = cbd_mix * frac_cbd;
        }
    }

    pub fn new(od: &OdMat, frac_self: FracSelf, frac_cbd:f64) -> MixMat {
        let mut mat = od.m.clone();
        let names = od.names.clone();
        let popns = od.popns.clone();
        Self::scale_diag(&mut mat, &frac_self);
        let cbd_ix = names.iter().position(|n| n == "20604").unwrap();
        Self::redistribute(&mut mat, &popns, frac_cbd, cbd_ix);
        MixMat {
            m: mat,
            names: names,
            popns: popns,
            frac_self: frac_self,
            frac_cbd: frac_cbd,
        }
    }
}

pub fn read_popns<P>(path: P, names: &Vec<String>) -> Array1<usize>
    where P: AsRef<Path> {
    let file = File::open(path).unwrap();
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b' ')
        .from_reader(file);
    let mut popn_tbl: HashMap<String, usize> = HashMap::new();
    for record in reader.records() {
        let row = record.unwrap();
        let name = row[1].to_string();
        if popn_tbl.contains_key(&name) {
            panic!("Multiple populations for {}", name);
        } else {
            let value = row[2].parse::<usize>().unwrap();
            popn_tbl.insert(name, value);
        }
    }
    let mut popns: Vec<usize> = vec![];
    for name in names.iter() {
        match popn_tbl.get(name) {
            Some(&popn) => popns.push(popn),
            None => panic!("No population defined for {}", name),
        }
    }
    Array::from_shape_vec(popns.len(), popns).unwrap()
}

pub fn find_mm<P>(path: P, subdir: bool) -> Vec<String>
    where P: AsRef<Path>{
    let mut files: Vec<String> = vec![];
    for entry in read_dir(path).unwrap() {
        let f = entry.unwrap().path();
        if f.is_file() {
            match f.extension() {
                Some(ext) if ext == "csv" =>
                    files.push(f.to_str().unwrap().into()),
                _ => ()
            }
        } else if f.is_dir() && subdir {
            for sub_entry in read_dir(f).unwrap() {
                let sub_f = sub_entry.unwrap().path();
                if sub_f.is_file() {
                    match sub_f.extension() {
                        Some(ext) if ext == "csv" =>
                            files.push(sub_f.to_str().unwrap().into()),
                        _ => ()
                    }
                }
            }
        }
    }
    files.sort();
    files
}
