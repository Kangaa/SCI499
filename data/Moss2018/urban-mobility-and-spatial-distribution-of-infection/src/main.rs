extern crate ndarray;
extern crate rayon;
extern crate serde;
extern crate spatial_seir;
#[macro_use]
extern crate serde_derive;
extern crate serde_json;
extern crate structopt;
#[macro_use]
extern crate structopt_derive;

use ndarray::prelude::*;
use rayon::prelude::*;
use spatial_seir::mixing::{FracSelf, OdMat, find_mm};
use spatial_seir::run::run_sims;
use std::fs::File;
use std::io::prelude::*;
use std::path::Path;
use structopt::StructOpt;

#[derive(Debug, Deserialize, Serialize)]
enum FracSelfVals {
    Uniform(Vec<f64>),
    Variable(Vec<Vec<f64>>),
}

#[derive(Debug, Deserialize, Serialize)]
struct Simulations {
    frac_self: FracSelfVals,
    frac_cbd: Vec<f64>,
    n_sims: usize,
    popn_file: String,
    mm_files: Vec<String>,
}

impl Default for Simulations {
    fn default() -> Simulations {
        Simulations {
            frac_self: FracSelfVals::Uniform(vec![0.1, 0.5, 0.75, 0.9]),
            frac_cbd: vec![0.20, 0.25, 0.33, 0.5],
            n_sims: 1000,
            popn_file: "./data/mixing/sa3-populations.ssv".to_string(),
            mm_files: find_mm("./data/mixing", true),
        }
    }
}

#[derive(Debug, StructOpt)]
#[structopt(name = "Spatial SEIR model", about = "Run epidemic simulations")]
struct Options {
    /// Whether to write each mixing matrix to disk.
    #[structopt(short = "w", long = "write-matrices",
                help = "Write mixing matrices to disk")]
    write_mats: bool,

    /// The optional simulation definition file.
    #[structopt(help = "Simulation definition file")]
    sim_file: Option<String>,
}

fn main () {
    let options = Options::from_args();

    let sims: Simulations = match options.sim_file {
        None => Default::default(),
        Some(ref sim_file) => {
            let sim_path = Path::new(&sim_file);
            if sim_path.exists() {
                let mut sim_file = File::open(&sim_path).unwrap();
                let mut sim_str = String::new();
                sim_file.read_to_string(&mut sim_str).unwrap();
                serde_json::from_str(&sim_str).unwrap()
            } else {
                println!("ERROR: file {:?} does not exist", sim_path);
                return;
            }
        },
    };

    let output_dir = "./data/model-output";
    let suffix = if sims.n_sims == 1 {
        "deterministic"
    } else {
        "output"
    };
    let frac_self = match sims.frac_self {
        FracSelfVals::Uniform(ref vals) => vals.iter()
            .map(|val| FracSelf::Uniform(*val))
            .collect(),
        FracSelfVals::Variable(ref vecs) => vecs.iter()
            .map(|vec| FracSelf::Variable(Array::from_vec(vec.clone())))
            .collect(),
    };
    let suffix = match sims.frac_self {
        FracSelfVals::Variable(_) => "variable-frac-self",
        _ => suffix,
    };

    sims.mm_files.par_iter()
        .for_each(|mm_file| {
            let od_mat = OdMat::read(mm_file, &sims.popn_file);
            let base = Path::new(&mm_file)
                .file_stem().unwrap()
                .to_str().unwrap();
            let out_path = format!("{}/{}-{}.ssv", output_dir, base, suffix);
            run_sims(&od_mat, &frac_self, &sims.frac_cbd, sims.n_sims,
                     &out_path, options.write_mats)});
}
