## TrialProbe Extraction Pipeline

This pipeline consists of Python code for processing ClinicalTrial.gov records into a simple JSON format.

In order to run this code, you need Python 3.9 and a set of dependencies specified in requirements.txt

Then, follow these steps:

1. Run the bash script download.sh to download all clinical trials with outcome data. `bash download.sh`
2. Run extract_raw.py to extract side effects and interventions and convert them to ICD10 and ATC codes. extract_raw.py requires one argument, a path to a download of a UMLS MRCONSO.RRF file, available from https://www.nlm.nih.gov/research/umls/licensedcontent/umlsknowledgesources.html. `python extract_raw.py --umls_mrconso_path ...`
3. Run extract_unique.py to combine trials that tested the same drugs and side effects and perform statistics calculations. `python extract_unique.py`

The result of this pipeline is a single file, unique_entries.txt. This file contains all the mappable clinical trial ground truth effects.

For the sake of convience, we have included a copy of that file in this repository that was computed from the March 3rd, 2021 ClinicalTrials.gov dump.
