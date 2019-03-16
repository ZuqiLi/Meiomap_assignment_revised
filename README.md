# Meiomap_assignment_revised
Revised code for Meiomap assignment

The revised code is structured in terms of a function, whose parameters are the input path of genotype file, the index of column used as reference and a boolean value which I'll explain later. By default, the index of reference column is set to be 6 which is the egg in the first trio. The function now is well generalised to multiple chromosomes and much more trios. Pandas and Numpy are applied in order to make the script efficient for large data and many trios.

The phasing method takes following strategies:

Homozygous maternal SNPs are non-informative but wi break regions of same phase (when appear within the region), labeled as -2.

NC SNPs in reference cell or any trio cells are also breakpoints of regions, labeled as -1.

The boolean parameter tells the function whether to show homozygous maternal SNPs and NC SNPs or not.

In the output file, two or more adjacent regions of the same phase are connected by breakpoints (shown or not).
