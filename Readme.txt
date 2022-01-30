Fortran Program: BP TAD Power Analysis.f
Function: to estimate the probability of n inversion breakpoints out of N total breakpoints being in a TAD.

Input files: 
boundary.tsv  - list of beginning and end coordinates of TAD boundaries for the third chromosome of D. pseudoobscura
Gene_Coor.tsv - list of beginning and end coordinates of genes on the third chromosome of D. pseudoobscura

Output files:
2021_10_10 BP_in_TADs_Bootstrap_Summary_Data n=13.txt 

Fortran Program: TAD Diff Expression Plot.f
Function: Draw a TAD map for the third chromosome of D. pseudoobscura along with inversion breakpoints

Input Files:
TADs.txt - list of beginning and end coordinates of TAD bodies for the third chromosome of D. pseudoobscura
BPs.txt  - list of beginning and end coordinates of inversion breakpoints for the third chromosome of D. pseudoobscura

Output Files:
TAD_Diff_Expression_Plot_Raw.ps - Adobe Postscript file of the TAD plot for Figure 3
Figure 3 TAD_Diff_Expression_Plot_Raw.pdf - Adobe pdf file created by conversion of the postscript file
