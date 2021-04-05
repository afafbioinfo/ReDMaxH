 #!/bin/sh

mkdir -p RedMaxHoutput
mkdir -p RedMaxPlots

# Make sure helix maximal generator is installed
# source https://github.gatech.edu/fhurley6/HelixEnumeration
 if [ ! -d ExternalScripts/HelixEnumeration-master/HelixEnumeration-master/Build ]; then
 

cd ExternalScripts/HelixEnumeration-master/HelixEnumeration-master
mkdir Build && cd Build
cmake ..
make
cd ..
cd ..
cd ..
cd ..
fi
pwd

echo "Parsing mutations...";
 #1- Parse shapemapper mut files
python2.7 ParseMapping.py --fasta Fasta/$1.fa  --Tag $2 OutputSHAPEMapper/Pipeline_Modified_$1_parsed.mut --untreated OutputSHAPEMapper/Pipeline_Untreated_$1_parsed.mut --mincoverage 0 $1
 
 echo "Mutations have been successfully parsed";

echo "Computing maximal helices...";
 
python2.7 GenerateMaximalHelices.py --Fasta Fasta/$1.fa --RNA $1 --ForwardAmpliconSize 0 --BackwardAmpliconSize 0 --minHelixsize 4
echo "Maximal helices have been successfully generated";

echo "Computing Relative differences on Moore neigborhood of order1";
python2.7 RedMaxH.py --file  $1FREQUENCY_$2ex.txt --output RedMaxHoutput --fasta Fasta/$1.fa   --a 0 --b  0 --Helixnbr 10 --p RedMaxHoutput/$1Profiling.csv

echo "Relative differences have been successfully computed";

echo "Generating plots...";
######## Generate heatmaps
Rscript --vanilla RelativeDifferenceViewer.R --RD1=./RedMaxHoutput/$1FREQUENCY_$2ex_Rawdata.csv --RD2=./RedMaxHoutput/$1FREQUENCY_$2ex_Rawdata.csv 
####### Clustering of ranked helices
Rscript --vanilla Clustering.R --C=./RedMaxHoutput/$1FREQUENCY_$2ex_contributingpairsforrankedhelices.csv  --R=./RedMaxHoutput/$1FREQUENCY_$2ex_RedmaxH_Rankedhelices.csv
echo "RedMaxH pipeline has been completed";
