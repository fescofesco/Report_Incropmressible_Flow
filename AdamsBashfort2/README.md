# Adams-Bashforth 2 solver

This directory contains the lecture-note-style AB2 incremental
pressure-correction implementation. It uses an Euler startup step, stores the
previous momentum and temperature RHS, includes the old pressure gradient in
the predictor, and performs one pressure projection per time step.

- Python entry point: `src_python/main.py`
- MATLAB entry point: `src_matlab/main.m`
- Generated results are written to this directory's `Plots_python` and
  `Plots_matlab` folders.

The folder name follows the requested repository layout. The standard method
name is "Adams-Bashforth" (with an `h`).
