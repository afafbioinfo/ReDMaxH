# RedMaxH
### mutation-Relative differences computed on Maximal Helices

 RedMaxH is a pipeline for predicting RNA structure features (=: Helices) compatible with DMS-MaP mutational data. From one or several input mutation data without any thermodynamic assumption, it computes and outputs the set of helices that are best supported by experimental data.

## Installing RedMaxH

RedMaxH consists in a set of Python 2.7+, C and R scripts. 


## Executing RedMaxH

RedMaxH can be invoked through a one bash command: 
      bash  5S DMSCellfree

      bash Processing.sh [ RNA] [condition]

**Example:** For an RNA sequence `5S`, and probing condition `DMSCellfree`, running the above command will lead to the production of outputs with the label `5S_DMSCellfree`.

The method will run with a default configuration.

## Input files


RedMaxH requires an RNA sequence with the SHAPEMapper files `.mut`.
 
 ### Sequence file
`{RNA}.fa` should be placed in the folder `Fasta`

### SHAPEMapper output file 
RedMaxH expects to find parsed files from SHAPEMapper `Pipeline_Untreated_/{RNA}/_parsed.mut` and `Pipeline_Modified_/{RNA}/_parsed.mut`, where `{RNA}` is the name of the chosen RNA (ie the name of the input FASTA file, minus its extension). 

Files should be placed in the folder `OutputSHAPEMapper`


## Outputs

###  Output Folders

RedMaxH typically produces many files in `RedMaxHoutput` and plots in `RedMaxPlots`

**RedMaxHoutput** files suffixe:

`Profiling.csv`  Set of maximal helices (name, opening, closing, length)

`Rawdata.csv`  for a tuple(i,j) [1-based], the MM count and the Rltdiff (=mutation relative difference) are reported. Note that we impose a minimal distance of 6 ncts 
between to count a joint mutation. 


`Norm_helices.csv`  Dataset ( type of the experiment: Cell-free, Incell...),Helix (Maximal helices),i (opening pair in the tuple (i,j) that belongs to the Moore neigborhood of the Helix),j (closing pair in the tuple (i,j) that belongs to the Moore neigborhood of the Helix),seqij (Nucleotide nature of the tuple(i,j),Rltdiff (Mutation relative difference),ijdist (length of the tuple)

`2DHelicesCommonPairs.csv`  Common tuples in the Moore neigborhood of order 1 between two helices.

`MergedHelices.csv `        List of helix pair compounds with the number of common tuples and the mutation relative difference of the coumpound.


`contributingpairsforrankedhelices.csv`  Value of relative difference for tuples from helix and helix coumpounds.

`RedmaxH_Rankedhelices.csv` The list of Ranked helices and helix-compounds 
(Helix (= Helix or coumpound) ,averageRLtdiff (Mean over all contributing tuples),nbrobservations (# of tuples),Length (Length of a Helix),nbrCCpairs,nbrAApairs,Profiling,std,median,averageonpositive,stdonpositive,medianonpsitive)


`Kmeans_elbow_mean.csv` Ranked helices affected to clusters ( where the number of clusters is determined using the elbow method) with their classification based on global std and global  mean.

**RedMaxPlots** files:


`Clusters_Helices.pdf` Clustering of Ranked helices according to  3 methods. If you are using the input example provided, you expect the result to be (link)![https://github.com/afafbioinfo/RedMaxH/blob/master/Clusters_Helices(Expectedplot).pdf]


