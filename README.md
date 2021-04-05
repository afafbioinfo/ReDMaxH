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

RedMaxHoutput
RedMaxPlots

IPANEMAP typically produces many messages during execution, to keep the user informed of its progress.
However, only the final (Pareto) structural models are output to the standard output, meaning that, after running

      python2.7 IPANEMAP.py > output.dat
      
the `output.dat` file will only consist of the final models.

**Example:** For an input sequence `GGGAAACCCAAAGGGAAACCC`, and probing profile assigning high accessibilities to `A`s, running the above command will lead to the production of a file `output.dat`, having content

    Structure                 dG   #SupportingConditions     BoltzmannProbability
    (((...)))...(((...)))   -4.3                       1       0.5735037115553154

where each line represents a cluster, and consists of:
  - Secondary structure model (centroid of the cluster)
  - Free-energy, as recomputed using `RNAeval`;
  - Number of supporting conditions;
  - Accumulated Boltzmann probability across conditions (aka stability in the companion manuscript), as computed using `RNAeval`. 
 
In this example, a unique probing condition implies a single model, but multiple structures may be produced in a multi-probing setting.

## Configuration
Most configuration options are set by modifying the content of the `IPANEMAP.cfg` file.

### Main options
 - `RNA`: Specifies a path (relative to the working directory) to a FASTA file where the nucleotide sequence of the main RNA of interest can be found. Note that the filename is important, as it will be used as a base name for the other input files. **Example:** `RNA: fasta_files/didymium.fa` will process the sequence found in the file, and `didymium` will be used as the *base name* of reactivities/hard contraints files (see `Conditions` option)
 - `SoftConstraintsDir` and `HardConstraintsDir`: Sets the **directories** used by IPANEMAP to locate soft (reactivities) and hard constraints files (if available)
 - `Conditions`: Can be used to specify the list of probing conditions used for the prediction. Should be set to a comma-separated list of conditions, i.e. the names of reactivity profiles/experiments to be considered for structure prediction
 
For an RNA having base name `{RNA}`, and a condition name `{Cond}`, IPANEMAP will attempt to locate files named `{SoftConstraintsDir}/{RNA}{Cond}.txt` and `{HardConstraintsDir}/{RNA}{Cond}.txt`. If none of these files is found, the method will rely on a purely thermodynamic sampling.

**Example:** Given a configuration
 
      [Input] 
      RNA: fasta_files/5sRNA.fa
      SoftConstraintsDir: soft
      HardConstraintsDir: hard
      Conditions: DMSMG,NMIA
      ...
   
the method will attempt to locate, and use for the sampling phase of the method, two files `5sRNADMSMG.txt` and `5sRNANMIA.txt` in each of the `soft` and `hard` directories.

### Sampling options
 - `DoSampling`: If set to `true`, IPANEMAP will always re-generate a representative structural sample (even if one can already be found)
 - `NumStructures`: Number of structures per condition, generated to approximate the pseudo-Boltzmann ensemble
 - `Temperature`: Temperature (in Celsius) used for the sampling
 - `m` and `b`: Slope and intercept used in the *reactivity to pseudo-energy* conversion (see Deigan et al, PNAS 2009)

### Misc options
 - `WorkingDir`: Main output directory for temp files, and final results of the analysis. Directory  will be created if non-existent.
 - `LogFile`: Name of file gathering the accumulated log. File will be created if non-existent.

### Visualization options
IPANEMAP currently relies on VARNA to produce
 - `DrawModels`: If set to `true`, uses VARNA to draw the final, Pareto-optimal, secondary structure models.
 - `DrawCentroids`: If set to `true`, uses VARNA to draw the centroids associated with all of the clusters.
 - `ShowProbing`:  If set to `true`, uses the reactivities of *the first probing condition* (as specified to the `cond` option, or  `Conditions` section of the config file) to annotate the secondary structure drawings.

Understand the fields of the ouput files:

Output Suffix:

Profiling.csv  Set of maximal helices (name, opening, closing, length)

Rawdata.csv  for a tuple(i,j) [1-based], the MM count and the Rltdiff (=mutation relative difference) are reported. Note that we impose a minimal distance of 6 ncts 
between to count a joint mutation. 


Norm_helices.csv  Dataset ( type of the experiment: Cell-free, Incell...),Helix (Maximal helices),i (opening pair in the tuple (i,j) that belongs to the Moore neigborhood of the Helix),j (closing pair in the tuple (i,j) that belongs to the Moore neigborhood of the Helix),seqij (Nucleotide nature of the tuple(i,j),Rltdiff (Mutation relative difference),ijdist (length of the tuple)

2DHelicesCommonPairs.csv  Common tuples in the Moore neigborhood of order 1 between two helices.

MergedHelices.csv         List of helix pair compounds with the number of common tuples and the mutation relative difference of the coumpound.


contributingpairsforrankedhelices.csv  Value of relative difference for tuples from helix and helix coumpounds.

RedmaxH_Rankedhelices.csv The list of Ranked helices and helix-compounds 
(Helix (= Helix or coumpound) ,averageRLtdiff (Mean over all contributing tuples),nbrobservations (# of tuples),Length (Length of a Helix),nbrCCpairs,nbrAApairs,Profiling,std,median,averageonpositive,stdonpositive,medianonpsitive)

