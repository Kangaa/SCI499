use ndarray::prelude::*;
use rand::{SeedableRng, StdRng, Rng};
use std::fs::File;
use std::path::Path;
use std::io::{Write, BufWriter};

use super::mixing::{FracSelf, MixMat, OdMat};
use super::seir::{Model, Params};

/// Create a new PRNG from a known seed.
pub fn rng_from_seed(seed: usize) -> StdRng {
    let seed_array: &[_] = &[seed];
    SeedableRng::from_seed(seed_array)
}

pub fn pick(v: Vec<f64>, rnd: f64) -> usize {
    let (x, _) = v.into_iter().enumerate().fold(
        (None, 0.0),
        |(m, sum), (ix, pr) | match m {
            Some(_) => (m, sum),
            None => {
                let sum = sum + pr;
                if sum > rnd {
                    (Some(ix), sum)
                } else {
                    (None, sum)
                }
            }
        },
    );
    x.unwrap()
}

pub struct RK4<'a> {
    k1: Array1<f64>,
    k2: Array1<f64>,
    k3: Array1<f64>,
    k4: Array1<f64>,
    k1_dt: Array1<f64>,
    k2_dt: Array1<f64>,
    k3_dt: Array1<f64>,
    x_k1: Array1<f64>,
    x_k2: Array1<f64>,
    x_k3: Array1<f64>,
    m: &'a Model<f64>,
    pub x: Array1<f64>,
    pub t: f64,
    dt: f64,
}

impl<'a> RK4<'a> {
    pub fn new(m: &Model<f64>, x0: Array1<f64>, t0: f64, dt: f64) -> RK4 {
        let shape = x0.len();
        RK4 {
            k1: Array::zeros(shape),
            k2: Array::zeros(shape),
            k3: Array::zeros(shape),
            k4: Array::zeros(shape),
            k1_dt: Array::zeros(shape),
            k2_dt: Array::zeros(shape),
            k3_dt: Array::zeros(shape),
            x_k1: Array::zeros(shape),
            x_k2: Array::zeros(shape),
            x_k3: Array::zeros(shape),
            m: m,
            x: x0,
            t: t0,
            dt: dt,
        }
    }

    pub fn step(&mut self) {
        let dt2 = self.dt / 2.0;
        
        self.k1 = self.m.f(&self.x);
        self.k1_dt = &self.k1 * dt2;
        self.x_k1 = &self.x + &self.k1_dt;

        self.k2 = self.m.f(&self.x_k1);
        self.k2_dt = &self.k2 * dt2;
        self.x_k2 = &self.x + &self.k2_dt;

        self.k3 = self.m.f(&self.x_k2);
        self.k3_dt = &self.k3 * self.dt;
        self.x_k3 = &self.x + &self.k3_dt;

        self.k4 = self.m.f(&self.x_k3);
        self.k2 *= 2.0;
        self.k3 *= 2.0;
        self.x += &((&self.k1 + &self.k2 + &self.k3 + &self.k4)
                    * self.dt / 6.0);
        self.t += self.dt;
    }
}

pub fn run_det<W>(params: Params, mat_info: &MixMat, rng: &mut StdRng,
                  stream: &mut BufWriter<W>)
    where W: Write {
    let infs = 10;
    let m = Model::<f64>::new(params, infs, rng);
    let t = 0.0 as f64;
    let dt = 0.1 as f64;

    // Run the ODE simulation to exhaustion.
    // Should instead consider a termination condition, such as having a net
    // prevalence (sum of all E and I compartments) less than one.
    let mut rk = RK4::new(&m, m.x0(), t, dt);
    while rk.t < 730.0 {
        rk.step();
    }

    // Report cumulative infections in each patch.
    let week = (rk.t / 7.0).floor() as usize;
    m.log_infs(&rk.x, &mat_info, week, stream);
}

pub fn run_stoch<W>(params: Params, mat_info: &MixMat, rng: &mut StdRng,
                  sim:usize, stream: &mut BufWriter<W>) -> bool
    where W: Write {
    let infs = 10;
    let mut m = Model::<usize>::new(params, infs, rng);
    let mut t = 0.0;
    let mut week = 0;

    loop {
        let rates = m.event_rates();
        let net_rate = rates.scalar_sum();
        if net_rate == 0.0 {
            break;
        }
        let dt = -rng.next_f64().ln() / net_rate;
        t += dt;
        // Note: suppress weekly output.
        if false && (t / 7.0).floor() as usize > week {
            // Report cumulative infections in each patch.
            m.log_infs(&mat_info, sim, week, stream);
            week += 1;
        }
        let rand_event = net_rate * rng.next_f64();
        let ix = pick(rates.to_vec(), rand_event);
        t += dt;
        m.event_occurred(ix);
    }

    if ! m.accept_simulation() {
        // Ignore this simulation.
        return false;
    }

    // Report cumulative infections in each patch.
    week = (t / 7.0).floor() as usize;
    m.log_infs(&mat_info, sim, week, stream);
    return true
}

pub fn run_sims<P>(od: &OdMat, f_self: &Vec<FracSelf>, f_cbd: &Vec<f64>,
                   n: usize, out_path: P, write_mats: bool)
    where P: AsRef<Path> {
    let mut rng = rng_from_seed(1234567);
    let r0 = 1.4;
    let sigma = 2.0;
    let gamma = 0.5;

    let path_ref = out_path.as_ref();
    println!("Writing to {} ...", path_ref.display());

    let mut stream = BufWriter::new(File::create(&path_ref).unwrap());
    writeln!(stream, "{} {} {} {} {} {} {} {}",
             "frac_self", "frac_cbd", "sim",
             "week", "ix", "cum_infs", "cum_inf_propn", "popn").unwrap();

    for ref frac_self in f_self.iter() {
        for &frac_cbd in f_cbd.iter() {
            let mat_info = MixMat::new(&od, (*frac_self).clone(), frac_cbd);
            let popns = &mat_info.popns;
            let mixmat = &mat_info.m;
            if write_mats {
                // Save this mixing matrix to disk.
                let frac_self_val = match frac_self {
                    &&FracSelf::Uniform(v) => v,
                    &&FracSelf::Variable(_) => 0.0,
                };
                let mat_file = format!(
                    "{}/{}-mixmat-{:.3}-{:.3}.ssv",
                    path_ref.parent().unwrap().to_str().unwrap(),
                    path_ref.file_stem().unwrap().to_str().unwrap(),
                    frac_self_val, frac_cbd);
                let mut m_stream = BufWriter::new(
                    File::create(&mat_file).unwrap());
                for row in mixmat.genrows() {
                    for value in row.iter() {
                        write!(m_stream, "{} ", value).unwrap();
                    }
                    writeln!(m_stream, "").unwrap();
                }
            }
            if n == 1 {
                let params = Params::new(r0, sigma, gamma, popns.clone(),
                                         mixmat.clone());
                run_det(params, &mat_info, &mut rng, &mut stream);
                continue
            }
            for sim in 0..n {
                let mut accepted = false;
                while ! accepted {
                    let params = Params::new(r0, sigma, gamma, popns.clone(),
                                             mixmat.clone());
                    accepted = run_stoch(params, &mat_info, &mut rng, sim,
                                         &mut stream);
                }
            }
        }
    }
}
