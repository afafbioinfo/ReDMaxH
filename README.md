# RedMaxH

ParseMapping allows to parse Shapemapper outputs and returns the raw mutation counts. It uses a cython code from Ringmapper program.

  - python ParseMapping.py --fasta RNA.fa Pipeline_Modified_RNA_parsed.mut --untreated  Pipeline_Untreated_RNA_parsed.mut   --mincoverage 0 OutputFile


Input files:

Insert Pipeline_Untreated_5S_parsed.mut and Pipeline_Modified_RNA_parsed.mut in the folder OutputSHAPEMapper 


Understand the ouput files:

Output Suffix:

Profiling.csv  Set of maximal helices (name, opening, closing, length)

Rawdata.csv  for a tuple(i,j) [1-based], the MM count and the Rltdiff (=mutation relative difference) are reported. Note that we impose a minimal distance of 6 ncts 
between to count a joint mutation. 


Norm_helices.csv  Dataset ( type of the experiment: Cell-free, Incell...),Helix (Maximal helices),i (opening pair in the tuple (i,j) that belongs to the Moore neigborhood of the Helix),j (closing pair in the tuple (i,j) that belongs to the Moore neigborhood of the Helix),seqij (Nucleotide nature of the tuple(i,j),Rltdiff (Mutation relative difference),ijdist (length of the tuple)

2DHelicesCommonPairs.csv  Common tuples in the Moore neigborhood of order 1 between two helices.

MergedHelices.csv         List of helix pair compounds with the number of common tuples and the mutation relative difference of the coumpound.


contributingpairsforrankedhelices.csv  Value of relative difference for tuples from helix and helix coumpounds.

RedmaxH_Rankedhelices.csv 
