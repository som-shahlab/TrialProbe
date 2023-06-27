# README

This code analyzes observational methods using the TrialProbe reference set.

plot.ipynb performs the core concordance analysis and plots both the concordance and recovery rate for entries when using the empirical Bayes filtering

compute_legacy_concordance.py performs our side analyasis that compares concordance with older approaches, including the OMOP, EU-ADR, and a variant of TrialProbe that is based on p values, not empirical Bayes.

The details of the conda environment necessary to run the notebooks and code is in environment.yaml.