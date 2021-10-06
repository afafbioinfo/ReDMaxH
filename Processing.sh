 #!/bin/sh

# Detect the système d exploitation
OS="`uname`"
case $OS in
  'Linux')
    OS='Linux'
    ##Generated shared library under linux
gcc -shared -o ExternalScripts/readMutStrings.so -fPIC ExternalScripts/readMutStrings.c

    ;;
  'Darwin') 
    OS='Mac'
    ##Generated shared library under OS X
g++ -dynamiclib -o ExternalScripts/readMutStrings.dylib ExternalScripts/readMutStrings.c

    ;;
  
  *) ;;
esac




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
python2.7 ParseMapping.py --fasta Fasta/$1.fa  --Tag $2 OutputSHAPEMapper/Pipeline_Modified_$1_parsed.mut  --mincoverage 0 $1 

echo "Mutations have been successfully parsed";

echo "Computing maximal helices...";
python2.7 GenerateMaximalHelices.py --Fasta Fasta/$1.fa --RNA $1 --ForwardAmpliconSize 0 --BackwardAmpliconSize 0 --minHelixsize 4

#144 mapped to 77
#python2.7 GenerateMaximalHelices.py --Fasta Fasta/$1.fa --RNA $1 --ForwardAmpliconSize 37 --BackwardAmpliconSize 31 --minHelixsize 4
#156 mapped to 77
#python2.7 GenerateMaximalHelices.py --Fasta Fasta/$1.fa --RNA $1 --ForwardAmpliconSize 49 --BackwardAmpliconSize 32 --minHelixsize 4

# 222 mapped to 77
#python2.7 GenerateMaximalHelices.py --Fasta Fasta/$1.fa --RNA $1 --ForwardAmpliconSize 115 --BackwardAmpliconSize 32 --minHelixsize 4

 
#156 mapped to 222 
#python2.7 GenerateMaximalHelices.py --Fasta Fasta/$1.fa --RNA $1 --Shift 207 --ForwardAmpliconSize 0 --BackwardAmpliconSize 0 --minHelixsize 4
#87 mapped to 222
# python2.7 GenerateMaximalHelices.py --Fasta Fasta/$1.fa --RNA $1 --Shift 264 --ForwardAmpliconSize 0 --BackwardAmpliconSize 0 --minHelixsize 4
#77 mapped to 222 
#python2.7 GenerateMaximalHelices.py --Fasta Fasta/$1.fa --RNA $1 --Shift 278 --ForwardAmpliconSize 0 --BackwardAmpliconSize 0 --minHelixsize 4
#144 mapped to 222 
#python2.7 GenerateMaximalHelices.py --Fasta Fasta/$1.fa --RNA $1 --Shift 217 --ForwardAmpliconSize 0 --BackwardAmpliconSize 0 --minHelixsize 4

#NYU 222 mapped to 156 python2.7 GenerateMaximalHelices.py --Fasta Fasta/$1.fa --RNA $1 --ForwardAmpliconSize 66 --BackwardAmpliconSize 0 --minHelixsize 4
## NYU222 python2.7 GenerateMaximalHelices.py --Fasta Fasta/$1.fa --RNA $1 --ForwardAmpliconSize 115 --BackwardAmpliconSize 32 --minHelixsize 4
## NYU156 python2.7 GenerateMaximalHelices.py --Fasta Fasta/$1.fa --RNA $1 --ForwardAmpliconSize 49 --BackwardAmpliconSize 32 --minHelixsize 4
### NUY144 python2.7 GenerateMaximalHelices.py --Fasta Fasta/$1.fa --RNA $1 --ForwardAmpliconSize 37 --BackwardAmpliconSize 31 --minHelixsize 4

###NUY87 python2.7 GenerateMaximalHelices.py --Fasta Fasta/$1.fa --RNA $1 --ForwardAmpliconSize 10 --BackwardAmpliconSize 0 --minHelixsize 4
echo "Maximal helices have been successfully generated";

echo "Computing Relative differences on Moore neigborhood of order1";
python2.7 RedMaxH.py --file  RedMaxHoutput/$1FREQUENCY_$2ex.txt --output RedMaxHoutput --fasta Fasta/$1.fa   --a 0 --b 0  --Helixnbr 10 --p RedMaxHoutput/$1Profiling.csv

echo "Relative differences have been successfully computed";

echo "Generating plots...";
######## Generate heatmaps
#Rscript --vanilla RelativeDifferenceViewer.R --RD1=./RedMaxHoutput/$1FREQUENCY_$2ex_Rawdata.csv --RD2=./RedMaxHoutput/$1FREQUENCY_$2ex_Rawdata.csv 
####### Clustering of ranked helices
Rscript --vanilla Clustering.R  --F $1 --C=./RedMaxHoutput/$1FREQUENCY_$2ex_contributingpairsforrankedhelices.csv  --R=./RedMaxHoutput/$1FREQUENCY_$2ex_RedmaxH_Rankedhelices.csv --P=./RedMaxHoutput/$1Profiling.csv --N=56
echo "RedMaxH pipeline has been completed";
