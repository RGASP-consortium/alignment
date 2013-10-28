This repository contains software for evaluation of spliced aligners,
written to assess RNA-seq mappers as part of the RGASP project.

Some of the scripts are tailored to particular data sets analyzed in
RGASP or to the computational environment at the European
Bioinformatics Institute, where the RGASP evaluation was carried out.

The evaluation considered alignments of both real and simulated
RNA-seq data. The simulated data was generated using the [BEERS]
simulator, and evaluated using the scripts in this repository. Some of
the scripts herein specifically deal with alignments of simulated
data, extracting accuracy metrics by comparison with the true
alignments and/or simulated transcript models. Other scripts extract
more general alignment statistics, and can be applied both to
alignments of real and simulated data. Most scripts that process
alignments expect input in SAM or BAM format.

A subset of the software written for the RGASP spliced aligner
evaluation are in a separate repository:
https://github.com/RGASP-consortium/coverage.
These scripts extract metrics relating to the coverage of genomic
features by alignments, and produce plots of those metrics.

Dependencies:
- [Genoman] - perl library required by some of the perl scripts
- [SAMTools] - for BAM file I/O
- R packages: gplots, plotrix, RColorBrewer

Please see further documentation under the doc directory.

If using this software, please cite:  
Engstr&ouml;m et al. (2013) Systematic evaluation of spliced alignment
programs for RNA-seq data. _Nat Methods_, in press.

If you have questions about this software, please write to
par.engstrom@scilifelab.se.

[BEERS]: http://cbil.upenn.edu/BEERS/
[Genoman]: http://www.ebi.ac.uk/~engstrom/genoman/
[SAMTools]: http://samtools.sourceforge.net/
