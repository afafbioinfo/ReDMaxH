
**Data**

4 RNAs:

5S rRNA E.Coli(Length =120, 80 nucleotides were considered (Length-2*Amplicon))

rpsB E.coli (Length=204, 164 ncts were considered)

rpsM E.coli (Lenght=241, 201 ncts were considered)

RMRP Homo Sapiens (Length=267, 227 ncts were considered)

**1/In Cell Vs.Cell Free**

We assume that the RNA folds the same way in Cell and in Cell-free condition. Is the difference between the distribution of the structure of tuples (i,j) statistically significant?



Let's start by performing KS test on 10.000 newly generated samples from the data(Permutation approach) and compare this KS-distribution to the kS-test of the data.


a- By combining the four RNAs, first we compare the **Cardinal of MM distribution between in cell and cell free conditions**:

| PP | UP | UU |
| :---         |     :---:      |          ---: |
| ![GitHub Logo](https://github.com/afafbioinfo/BackupGatech/blob/master/April2nd2020/1-Figures/PP_Cell_CellFree_4RNAs.png)|![GitHub Logo](https://github.com/afafbioinfo/BackupGatech/blob/master/April2nd2020/1-Figures/UP_Cell_CellFree_4RNAs.png)|![GitHub Logo](https://github.com/afafbioinfo/BackupGatech/blob/master/April2nd2020/1-Figures/UU_Cell_CellFree_4RNAs.png)   |
| 16544 Obs. D = 0.021035, p-value = 0.001324    | 28418 Obs. D = 0.044866, p-value < 2.2e-16       | 11344 Obs. D = 0.068671, p-value < 2.2e-16   |

b- Then we compare the Cardinal MM and the Prior accross the two conditions:
|PP Cell-Free | PP IN cell |
| :---         |     :---:      |    
|![GitHub Logo](https://github.com/afafbioinfo/BackupGatech/blob/master/April2nd2020/1-Figures/CardinalvsPriorPP_CellFree_4RNAs.png)|![GitHub Logo](https://github.com/afafbioinfo/BackupGatech/blob/master/April2nd2020/1-Figures/CardinalvsPriorPP_Cell_4RNAs.png)|
|D = 0.024178, p-value = 0.0001262 | D = 0.030404, p-value = 4.564e-07|


|UP Cell-Free  | UP IN cell|
| :---         |     :---:      |    
|![GitHub Logo](https://github.com/afafbioinfo/BackupGatech/blob/master/April2nd2020/1-Figures/CardinalvsPriorUP_CellFree_4RNAs.png)|![GitHub Logo](https://github.com/afafbioinfo/BackupGatech/blob/master/April2nd2020/1-Figures/CardinalvsPriorUP_Cell_4RNAs.png)|
|D = 0.0090436, p-value = 0.1955 | D = 0.0099233, p-value = 0.1218|



|UU Cell-Free | UU IN Cell|
| :---         |     :---:      |    
|![GitHub Logo](https://github.com/afafbioinfo/BackupGatech/blob/master/April2nd2020/1-Figures/CardinalvsPriorUU_CellFree_4RNAs.png)|![GitHub Logo](https://github.com/afafbioinfo/BackupGatech/blob/master/April2nd2020/1-Figures/CardinalvsPriorUU_Cell_4RNAs.png)|
|D = 0.0047602, p-value = 0.9995 | D = 0.0076693, p-value = 0.8925|

*Remarks*

*When comparing data for 4 RNAs from Cell and Cell-free conditions, the Prior seems to be a good approximation for (i,j) in [(Upaired,Unpaired),(Unpaired,paired),(paired,Unpaired)] with a P-value above 0.05.

** 2/A nucleotide and a structure decomposition of the mutation rate distributions **

DMS probes the base, what does the Mutation rate distribution looks like for each single nucleotides(C,G,U,A) and for each structural context (Unpaired, Helix-end, paired in Stack)?

The obtained distributions for both Cell and Cell-Free conditions sorted by nucleotide nature:

|Cell-Free | IN Cell|
| :---         |     :---:      |    
|![GitHub Logo](https://github.com/afafbioinfo/BackupGatech/blob/master/April2nd2020/2--Figures/MutationrateC-A-CellFree.png)|![GitHub Logo](https://github.com/afafbioinfo/BackupGatech/blob/master/April2nd2020/2--Figures/MutationrateC-A-INCell.png)|
|![GitHub Logo](https://github.com/afafbioinfo/BackupGatech/blob/master/April2nd2020/2--Figures/MutationrateG-U-CellFree.png)|![GitHub Logo](https://github.com/afafbioinfo/BackupGatech/blob/master/April2nd2020/2--Figures/MutationrateG-U-INCell.png)|


The obtained distributions for both Cell and Cell-Free conditions sorted by local structure:
    
298 Obs. for Unpaired, 242 Obs. for Stack paired and 121 Obs. for Helix-end paired nucleotides.

|Cell-Free | IN Cell|
| :---         |     :---:      |
|![GitHub Logo](https://github.com/afafbioinfo/BackupGatech/blob/master/April2nd2020/2--Figures/Unpaired_HE_Stack_CELLFree.png)|![GitHub Logo](https://github.com/afafbioinfo/BackupGatech/blob/master/April2nd2020/2--Figures/Unpaired_HE_Stack_INCELL.png)|


Let us decompose these distributions in function of the nucleotide nature:
|Nucleotide|Cell-Free | IN Cell|
| :---    | :---         |     :---:      |
|A|![GitHub Logo](https://github.com/afafbioinfo/BackupGatech/blob/master/April2nd2020/2--Figures/A-Unpaired_HE_Stack_CELLFree.png)|![GitHub Logo](https://github.com/afafbioinfo/BackupGatech/blob/master/April2nd2020/2--Figures/A-Unpaired_HE_Stack_INCELL.png)|
|C|![GitHub Logo](https://github.com/afafbioinfo/BackupGatech/blob/master/April2nd2020/2--Figures/C-Unpaired_HE_Stack_CELLFree.png)|![GitHub Logo](https://github.com/afafbioinfo/BackupGatech/blob/master/April2nd2020/2--Figures/C-Unpaired_HE_Stack_INCELL.png)|
|G|![GitHub Logo](https://github.com/afafbioinfo/BackupGatech/blob/master/April2nd2020/2--Figures/G-Unpaired_HE_Stack_CELLFree.png)|![GitHub Logo](https://github.com/afafbioinfo/BackupGatech/blob/master/April2nd2020/2--Figures/G-Unpaired_HE_Stack_INCELL.png)|
|U|![GitHub Logo](https://github.com/afafbioinfo/BackupGatech/blob/master/April2nd2020/2--Figures/U-Unpaired_HE_Stack_CELLFree.png)|![GitHub Logo](https://github.com/afafbioinfo/BackupGatech/blob/master/April2nd2020/2--Figures/U-Unpaired_HE_Stack_INCELL.png)|

3/ 5SrRNA_Ecoli:
 
A decomposition of the distribution of the difference based on the nucleotide nature:
| :---         |     :---:      |
|(i,j) has value in [A]|![GitHub Logo](https://github.com/afafbioinfo/BackupGatech/blob/master/April2nd2020/3--Figures/AA/5S_cellfree_5S_DMSmodified.txt_%20%23MM(i%2Cj)-Mutr(i)*Mutr(j)*Nij.png)|
|(i,j) has value in [C]|[GitHub Logo](https://github.com/afafbioinfo/BackupGatech/blob/master/April2nd2020/3--Figures/CC/5S_cellfree_5S_DMSmodified.txt_%20%23MM(i%2Cj)-Mutr(i)*Mutr(j)*Nij.png)|
|(i,j) has value in [G]|[GitHub Logo](https://github.com/afafbioinfo/BackupGatech/blob/master/April2nd2020/3--Figures/GG/5S_cellfree_5S_DMSmodified.txt_%20%23MM(i%2Cj)-Mutr(i)*Mutr(j)*Nij.png)|
|(i,j) has value in [U]|[GitHub Logo](https://github.com/afafbioinfo/BackupGatech/blob/master/April2nd2020/3--Figures/UU/5S_cellfree_5S_DMSmodified.txt_%20%23MM(i%2Cj)-Mutr(i)*Mutr(j)*Nij.png)|




(i,j) has value in [C,A]
3--Figures/CA/5S_cellfree_5S_DMSmodified.txt_ #MM(i,j)-Mutr(i)*Mutr(j)*Nij.png
(i,j) has value in [U,G]
3--Figures/UG/5S_cellfree_5S_DMSmodified.txt_ #MM(i,j)-Mutr(i)*Mutr(j)*Nij.png
(i,j) has value in [G,C]
3--Figures/GC/5S_cellfree_5S_DMSmodified.txt_ #MM(i,j)-Mutr(i)*Mutr(j)*Nij.png
(i,j) has value in [U,A]
3--Figures/UA/5S_cellfree_5S_DMSmodified.txt_ #MM(i,j)-Mutr(i)*Mutr(j)*Nij.png
(i,j) has value in [G,A]
3--Figures/GA/5S_cellfree_5S_DMSmodified.txt_ #MM(i,j)-Mutr(i)*Mutr(j)*Nij.png
(i,j) has value in [C,U]
3--Figures/CU/5S_cellfree_5S_DMSmodified.txt_ #MM(i,j)-Mutr(i)*Mutr(j)*Nij.png




\bibliographystyle{plain}
%\bibliography{references}
\end{document}
