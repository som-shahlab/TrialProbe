## TrialProbe

**TrialProbe** is a toolkit for probing non-experimental methods for inferring cauasal effects using ClinicalTrials.gov data.

This codebase consists of two components which correspond to two folders, extraction and analysis.

*create_reference_set* is Python and Julia code for extracting, cleaning, and preprocessing ClinicalTrials.gov data into a simple JSON format, reference_set.txt

To encourage ease of use, an example reference_set.txt from that pipeline can be found in this folder.

*analyze_methods* is Python code for analyzing non-experimental methods using reference_set.txt from the *create_refererence_set* pipeline.

The recommended entry point for new users is analyze_methods/plot.ipynb, as it illustrates how to use this reference set to compute concordance.