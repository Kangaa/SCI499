//! # SEIR model
//!
//! ## Overview
//!
//! ```rust
//! # extern crate rand;
//! # extern crate spatial_seir;
//! #
//! # fn main() {
//! use rand::Rng;
//! use std::fs::File;
//! use std::path::Path;
//! use std::io::{Write, BufWriter};
//! use spatial_seir::seir;
//! use spatial_seir::mixing;
//! use spatial_seir::run;
//!
//! // Population mixing characteristics and patch populations.
//! let frac_self = 0.5;
//! let frac_cbd = 0.2;
//! let path = "data/mixing/pub_PR_Normed.csv";
//! let mat_info = mixing::read_mat(&path, frac_self, frac_cbd);
//! let mixmat = mat_info.m.clone();
//! let names = mat_info.names.clone();
//! let popns = mat_info.popns.clone();
//!
//! // Epidemiological parameters.
//! let r0 = 1.4;
//! let sigma = 2.0;
//! let gamma = 0.5;
//!
//! // Model parameters.
//! let params = seir::Params::new(r0, sigma, gamma, popns, mixmat);
//!
//! // Create the output file. It will contain the following columns:
//! // - frac_self: the value of the `frac_self` parameter.
//! // - frac_cbd: the value of the `frac_cbd` parameter.
//! // - sim: an integer, used to uniquely identify each simulation.
//! // - week: the week number (0..n).
//! // - patch: the patch name/identifier.
//! // - cum_infs: the cumulative number of infections in the patch.
//! // - cum_inf_propn: the proportion of the population in the patch that
//! //   have been infected.
//! // - popn: the patch population.
//! let mut output = BufWriter::new(File::create("output.ssv").unwrap());
//! writeln!(output, "{} {} {} {} {} {} {} {}",
//!          "frac_self", "frac_cbd", "sim", "week", "patch",
//!          "cum_infs", "cum_inf_propn", "popn").unwrap();
//!
//! // Simulate model events until the epidemic is finished.
//! let mut rng = spatial_seir::run::rng_from_seed(1001);
//! let sim = 0;
//! run::run(params, &mat_info, &mut rng, sim, &mut output);
//! # }
//! ```

use ndarray::prelude::*;
use rand::{StdRng, Rng};
use std::io::{Write, BufWriter};

use super::mixing::{FracSelf, MixMat};

/// Model parameters.
pub struct Params {
    /// The daily force of infection.
    pub beta: f64,
    /// The infectiousness rate.
    pub sigma: f64,
    /// The recovery rate.
    pub gamma: f64,
    /// The resident population in each patch.
    pub popn: Array1<usize>,
    /// The mixing matrix.
    pub mixmat: Array2<f64>,
}

impl Params {
    pub fn new(r0: f64, sigma: f64, gamma: f64, popn: Array1<usize>,
               mixmat: Array2<f64>) -> Params {
        if ! mixmat.is_square() {
            panic!("Mixing matrix is not square")
        }
        if mixmat.cols() != popn.len() {
            panic!("Population vector does not match mixing matrix")
        }
        Params {
            beta: r0 * gamma,
            sigma: sigma,
            gamma: gamma,
            popn: popn,
            mixmat: mixmat,
        }
    }
}

/// Model state.
pub struct Model<T> {
    /// The model parameters.
    pub params: Params,
    /// The number of patches.
    pub np: usize,
    /// The current time.
    pub t: f64,
    /// The susceptible population in each patch.
    pub s: Array1<T>,
    /// The exposed population in each patch.
    pub e: Array1<T>,
    /// The infectious population in each patch.
    pub i: Array1<T>,
}

/// Model events.
pub enum Event {
    /// An exposure event in patch `i`.
    ///
    /// ```rust
    /// use spatial_seir::seir::Event;
    /// let i = 0;
    /// let ev = Event::Exposure(i);
    /// ```
    Exposure(usize),
    /// An infection event in patch `i`.
    ///
    /// ```rust
    /// use spatial_seir::seir::Event;
    /// let i = 6;
    /// let ev = Event::Infection(i);
    /// ```
    Infection(usize),
    /// A recovery event in patch `i`.
    ///
    /// ```rust
    /// let i = 2;
    /// use spatial_seir::seir::Event;
    /// let ev = Event::Recovery(i);
    /// ```
    Recovery(usize),
}

impl Model<f64> {
    pub fn new(params: Params, i0: usize, rng: &mut StdRng) -> Model<f64> {
        let n_patch = params.popn.len();
        let rand_choice = rng.next_f64();
        let ix = (rand_choice * (n_patch as f64)).trunc() as usize;
        assert!(ix < n_patch, "Selected patch {} of {}, from {}", ix, n_patch,
                rand_choice);
        // let s0 = params.popn.to_owned();
        let s0 = params.popn.mapv(|x| x as f64);
        let mut m: Model<f64> = Model {
            params: params,
            np: n_patch,
            t: 0.0,
            s: s0,
            e: Array::zeros(n_patch),
            i: Array::zeros(n_patch),
        };
        m.s[ix] -= i0 as f64;
        m.e[ix] += i0 as f64;
        m
    }

    pub fn log_infs<W>(&self, x: &Array1<f64>, mat_info: &MixMat, week: usize,
                       stream: &mut BufWriter<W>) where W: Write {
        let n1 = self.np as isize;
        let n2 = 2 * self.np as isize;
        let s = x.slice(s![..n1]);
        let e = x.slice(s![n1..n2]);
        let cum_infs = self.params.popn.mapv(|x| x as f64) - (&s + &e);
        assert!(cum_infs.len() == mat_info.names.len(),
                "Have {} patches and {} patch names",
                cum_infs.len(), mat_info.names.len());
        assert!(cum_infs.len() == self.params.popn.len(),
                "Have {} patches and {} populations",
                cum_infs.len(), self.params.popn.len());
        let frac_self = match mat_info.frac_self {
            FracSelf::Uniform(val) => val,
            FracSelf::Variable(_) => 0.0,
        };
        for (ix, cinf) in cum_infs.iter().enumerate() {
            let propn = cinf / self.params.popn[ix] as f64;
            let name = &mat_info.names[ix];
            let popn = self.params.popn[ix];
            let sim = 0;
            writeln!(stream, "{} {} {} {} {} {} {} {}",
                     frac_self, mat_info.frac_cbd, sim, week, name,
                     cinf, propn, popn).unwrap();
            stream.flush().unwrap();
        }
    }

    pub fn rates(&self) -> Array1<f64> {
        let mut rates = Array::zeros(3 * self.params.popn.len());
        let s_frac: Array1<f64> = &self.s
            / &self.params.popn.map(|&x| x as f64);
        // Calculate the force of infection, allowing for mixing.
        let inf_force = self.params.beta * self.params.mixmat.t().dot(&self.i);
        let n1 = self.np as isize;
        let n2 = 2 * self.np as isize;
        let exp_rate = s_frac * inf_force;
        let inf_rate = self.params.sigma * &self.e;
        let rec_rate = self.params.gamma * &self.i;
        let d_s = -&exp_rate;
        let d_e = &exp_rate - &inf_rate;
        let d_i = &inf_rate - &rec_rate;
        rates.slice_mut(s![..n1]).assign(&d_s);
        rates.slice_mut(s![n1..n2]).assign(&d_e);
        rates.slice_mut(s![n2..]).assign(&d_i);
        rates
    }

    pub fn x0(&self) -> Array1<f64> {
        let mut x0 = Array::zeros(3 * self.params.popn.len());
        let n1 = self.np as isize;
        let n2 = 2 * self.np as isize;
        x0.slice_mut(s![..n1]).assign(&self.s);
        x0.slice_mut(s![n1..n2]).assign(&self.e);
        x0.slice_mut(s![n2..]).assign(&self.i);
        x0
    }

    pub fn f(&self, x: &Array1<f64>) -> Array1<f64> {
        let exp_len = 3 * self.params.popn.len();
        let n1 = self.np as isize;
        let n2 = 2 * self.np as isize;
        let mut rates = Array::zeros(exp_len);
        assert!(x.len() == exp_len, "Expected length {}, and {}",
                rates.len(), x.len());
        let s = x.slice(s![..n1]);
        let e = x.slice(s![n1..n2]);
        let i = x.slice(s![n2..]);
        let s_frac: Array1<f64> = &s
            / &self.params.popn.map(|&x| x as f64);
        // Calculate the force of infection, allowing for mixing.
        let inf_force = self.params.beta * self.params.mixmat.t().dot(&i);
        let n1 = self.np as isize;
        let n2 = 2 * self.np as isize;
        let exp_rate = s_frac * inf_force;
        let inf_rate = self.params.sigma * &e;
        let rec_rate = self.params.gamma * &i;
        let d_s = -&exp_rate;
        let d_e = &exp_rate - &inf_rate;
        let d_i = &inf_rate - &rec_rate;
        rates.slice_mut(s![..n1]).assign(&d_s);
        rates.slice_mut(s![n1..n2]).assign(&d_e);
        rates.slice_mut(s![n2..]).assign(&d_i);
        rates
    }
}

impl Model<usize> {
    pub fn new(params: Params, i0: usize, rng: &mut StdRng) -> Model<usize> {
        let n_patch = params.popn.len();
        let rand_choice = rng.next_f64();
        let ix = (rand_choice * (n_patch as f64)).trunc() as usize;
        assert!(ix < n_patch, "Selected patch {} of {}, from {}", ix, n_patch,
                rand_choice);
        let s0 = params.popn.to_owned();
        let mut m = Model {
            params: params,
            np: n_patch,
            t: 0.0,
            s: s0,
            e: Array::zeros(n_patch),
            i: Array::zeros(n_patch),
        };
        m.s[ix] -= i0;
        m.e[ix] += i0;
        m
    }

    pub fn event_rates(&self) -> Array1<f64> {
        let mut rates = Array::zeros(3 * self.params.popn.len());
        let s_frac: Array1<f64> = &self.s.map(|&x| x as f64)
            / &self.params.popn.map(|&x| x as f64);
        let infs: Array1<f64> = self.i.map(|&x| x as f64);
        // Calculate the force of infection, allowing for mixing.
        let inf_force = self.params.beta * self.params.mixmat.t().dot(&infs);
        let n1 = self.np as isize;
        let n2 = 2 * self.np as isize;
        let exp_rate = s_frac * inf_force;
        rates.slice_mut(s![..n1]).assign(&exp_rate);
        let inf_rate = self.params.sigma * &self.e.map(|&x| x as f64);
        rates.slice_mut(s![n1..n2]).assign(&inf_rate);
        let rec_rate = self.params.gamma * self.i.map(|&x| x as f64);
        rates.slice_mut(s![n2..]).assign(&rec_rate);
        rates
    }

    pub fn event_occurred(&mut self, event_ix: usize) {
        let ev_type: usize = event_ix / self.np;
        let ev_locn: usize = event_ix % self.np;
        if ev_type == 0 {
            self.s[ev_locn] -= 1;
            self.e[ev_locn] += 1;
        } else if ev_type == 1 {
            self.e[ev_locn] -= 1;
            self.i[ev_locn] += 1;
        } else if ev_type == 2 {
            self.i[ev_locn] -= 1;
        } else {
            panic!("Invalid event type {}", ev_type);
        }
    }

    pub fn log_infs<W>(&self, mat_info: &MixMat, sim: usize, week: usize,
                       stream: &mut BufWriter<W>) where W: Write {
        let cum_infs = self.params.popn.mapv(|x| x as f64) -
            (self.s.map(|&x| x as f64) + self.e.map(|&x| x as f64));
        assert!(cum_infs.len() == mat_info.names.len(),
                "Have {} patches and {} patch names",
                cum_infs.len(), mat_info.names.len());
        assert!(cum_infs.len() == self.params.popn.len(),
                "Have {} patches and {} populations",
                cum_infs.len(), self.params.popn.len());
        let frac_self = match mat_info.frac_self {
            FracSelf::Uniform(val) => val,
            FracSelf::Variable(_) => 0.0,
        };
        for (ix, cinf) in cum_infs.iter().enumerate() {
            let propn = cinf / self.params.popn[ix] as f64;
            let name = &mat_info.names[ix];
            let popn = self.params.popn[ix];
            writeln!(stream, "{} {} {} {} {} {} {} {}",
                     frac_self, mat_info.frac_cbd, sim, week, name,
                     cinf, propn, popn).unwrap();
            stream.flush().unwrap();
        }
    }

    /// Check whether the epidemic was sufficiently large.
    pub fn accept_simulation(&self) -> bool {
        let cum_infs = self.params.popn.mapv(|x| x) -
            (self.s.map(|&x| x) + self.e.map(|&x| x));
        cum_infs.scalar_sum() >= 50000
    }
}
