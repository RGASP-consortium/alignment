This repository contains software for evaluation of spliced aligners,
written to assess RNA-seq mappers as part of the RGASP project.

Some of the scripts are tailored to particular data sets analyzed in
RGASP or to the computational environment at the European
Bioinformatics Institute, where the RGASP evaluation was carried out.
We are working to generalize the code and expand the documentation to
simplify reuse for other projects.

The evaluation considered alignments of both real and simulated
RNA-seq data. The simulated data was generated using the [BEERS]
simulator, and evaluated using the scripts in this repository. Some of
the scripts herein specifically deal with alignments of simulated
data, extracting accuracy metrics by comparison with the true
alignments and/or simulated transcript models. Other scripts extract
more general alignments statistics, and can be applied both to
alignments of real and simulated data. Most scripts that process
alignments expect input in SAM or BAM format.

A subset of the software written for the RGASP spliced aligner
evaluation are in a separate repository:
https://github.com/sbotond/paper-rgasp3-cov.
These scripts extract metrics relating to the coverage of genomic
features by alignments, and produce plots of those metrics.

Please see further documentation under the doc directory.

If using this software, please cite:  
Engstr&ouml;m et al. Systematic evaluation of spliced aligners for RNA-seq, _submitted_.

[BEERS]: http://cbil.upenn.edu/BEERS/
