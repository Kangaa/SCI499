# Overview

This directory contains the simulation model and analysis scripts to produce
all of the manuscript figures.

# Usage

Analyse the simulation outputs in the `data/model-output` directory and plot
the results:

    make

# Requirements

- R 3.4 and the following packages:

  - Cairo
  - ggplot2
  - gridExtra
  - RColorBrewer
  - viridis

# Simulation model

The simulation model requires:

- [Rust 1.21](https://blog.rust-lang.org/2017/10/12/Rust-1.21.html) or later.
  The preferred method to install Rust is via [rustup](https://www.rustup.rs/).

To run the model simulations and overwrite the existing output files in the
`data/model-output` directory:

    cargo run --release -- sim/default-ode.json
