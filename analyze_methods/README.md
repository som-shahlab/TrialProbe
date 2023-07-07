# TrialProbe Method Analysis Pipeline

This code analyzes non-experimental methods using the TrialProbe reference set.

# Setup

This code requires Python 3.10
We have provided an environment.yml conda environment with the required dependencies.

conda env create -n trialprobe --file environment.yml

# Steps to Run

You must first run non-experimental methods on each item in the reference set (see ../reference_set.txt) and store the results in the required format demonstrated in observational_results.txt

observational_results.txt contains several non-experimental methods computed on Optum data.

plot.ipynb then performs the core concordance analysis and plots both the concordance and recovery rate for entries at various effect size thresholds.

compute_legacy_concordance.py performs our side analyasis that compares concordance with older approaches, including the OMOP, EU-ADR, and a variant of TrialProbe that is based on Fisher exact tests.

create_table.py, create_table2.py, create_smaller_conso.py, plot_reference_set_stats.ipynb and compare_reference_sets.ipynb are the code to reproduce figures and analysis from our paper, but are probably not of general use to other users.
