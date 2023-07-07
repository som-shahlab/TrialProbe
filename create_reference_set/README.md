## TrialProbe Create Reference Set Pipeline

This pipeline consists of Python code for processing ClinicalTrial.gov records into a simple JSON format.
It also performs empirical Bayes analysis for estimating the effects in the trial records.

# Setup

This code requires Python 3.10 and Julia. 
We have provided an environment.yml conda environment with the required dependencies.

conda env create -n trialprobe_extraction --file environment.yml

This will install Python, Julia and all the necessary Python packages.

You must also install the Julia packages spereately with the following commands in the julia command line

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

# Steps to Run

Running this pipeline reuires the following steps

```
bash download.sh
python extract_raw.py --umls_mrconso_path ...
python extract_unique.py --require_high_quality_rcts
julia compute_statistics.jl
```

The output of this program is a single file, reference_set.txt with all the extracted effects and statistical analysis.

We have included a reference copy of reference_set.txt using a June 4rd, 2020 download of ClinicalTrials.gov in the main TrialProbe folder.
