//! # Spatially explicit SEIR model of influenza transmission.
//!

extern crate csv;
#[macro_use(s)]
extern crate ndarray;
extern crate rand;
extern crate rayon;

pub mod mixing;
pub mod seir;
pub mod run;
