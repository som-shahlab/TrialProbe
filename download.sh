wget https://clinicaltrials.gov/AllPublicXML.zip
mkdir nct
unzip AllPublicXML.zip -d nct

LC_ALL=C grep -F "<reported_events>" -r nct > nct_with_events.txt
