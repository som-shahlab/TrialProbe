# README

## R code

The file `json_to_csv.R` takes as input the JSON file `results.json` and converts it to the CSV file `trialverify_full.csv`.

## Julia code

The file `deconvolve.jl` conducts the deconvolution, the evaluation, and generates all plots.

The above code was run on Julia version 1.7.2. All required Julia packages and their version are specified in `Project.toml` and `Manifest.toml`. They may be installed automatically by starting a Julia session in this folder and typing:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```