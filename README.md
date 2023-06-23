## TrialProbe

**TrialProbe** is a toolkit for probing non-experimental methods for inferring cauasal effects using ClinicalTrials.gov data.

This codebase consists of two components which correspond to two folders, extraction and analysis.

*extraction* is Python code for extracting, cleaning, and preprocessing ClinicalTrials.gov data into a simple JSON format

*analysis* is Julia code for analyzing non-experimental methods using the proccessed JSON data from the *extraction* pipeline

To encourage ease of use, we have also stored the output of the *extraction* pipeline within the *analysis* folder so users can avoid reprocessing ClinicalTrials.gov in order to analyze a new method.

