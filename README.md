## TrialProbe

**TrialProbe** is a toolkit for probing non-experimental methods for inferring cauasal effects using ClinicalTrials.gov data.

This codebase consists of two components which correspond to two folders, extraction and analysis.

*extraction* is Python and Julia code for extracting, cleaning, and preprocessing ClinicalTrials.gov data into a simple JSON format, reference_set.txt

*analysis* is Python code for analyzing non-experimental methods using reference_set.txt from the *extraction* pipeline.

To encourage ease of use, we have also stored the output of the *extraction* pipeline within this folder as reference_set.txt so users can avoid reprocessing ClinicalTrials.gov in order to analyze a new method.

The recommended entry point for new users is analyis/plot.ipynb, as it illustrates how to use this reference set to compute concordance.