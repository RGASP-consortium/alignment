Manual: RGASP alignment evaluation scripts
==========================================

This text describes how to use the alignment evaluation scripts
written for RGASP, starting from BAM files to obtain tables and plots
of results.

Note: this documentation is a work in progress and not yet complete!
In the mean time, please send any questions to par.engstrom@scilifelab.se.

Prerequisites
-------------

First, make sure the dependencies are installed (see README.md). Then
set following two environment variable, preferably in your
~/.bash_profile (or the equivalent file if you are using a shell other than bash):

- RGASP_ALI_HOME  Base directory of this package
- RGASP_ALI_DATA  Data directory to use (can be your own data or the sample data included)

For example:

    RGASP_ALI_HOME=~/rgasp/alignment
    RGASP_ALI_DATA=~/rgasp/alignment/data

Organization of scripts
-----------------------

The processing of alignment files to extract statistics is carried out
by a number of perl scripts. The repository also includes several
R scripts that read the output from the perl scripts, and summarize
the data in plots and tables. (Note: plotting scripts have not yet
been added to repository.)

Most of the perl scripts process a single alignment file, and thus do
not make assumptions about the naming of files or their organization
in the filesystem. Exceptions are validate_bam.pl and
cufflinks_eval.pl (but we shall strive to generalize those further).

In contrast, the R scripts generally summarize results collected from
different alignment files to present them for comparison. These
scripts therefore use configuration files that list aligners and data
sets, and make assumptions about where data files are located. All
configuration and data files should be in subdirectories under the
path specified in $RGASP_ALI_DATA.

Configuration files
-------------------

The R scripts use three configuration files:

- $RGASP_ALI_DATA/aligners/methods.txt
- $RGASP_ALI_DATA/aligners/styles.txt
- $RGASP_ALI_DATA/reads/datasets.txt

### methods.txt ###

This file lists the evaluated alignment protocols.

- name: Symbolic name for protocol. May contain spaces.
- team: Name of developer team contributing the protocol. Should
  consist of word characters only (A-z, 0-9 and _).
- number: Protocol number (a team may contribute several protocols and each should have a unique number)
- description: Description of protocol (optional)
- color: color used for plotting symbols associated with this protocol
- style: plotting style used for this protocol (see styles.txt)

Note that the combination of _team_ and _number_ is used as a unique ID for
each protocol, and _name_ should also be unique.

The repository currently contains a file methods_all.txt with all
protocols evaluated in RGASP, as well as a file methods.txt with a
subset of protocols, for which data files are provided. (Note: data
files have not yet been added to the repository. Following acceptance
of the manuscript we shall add data files for all protocols, and prior
to that for the subset listed in methods.txt.)

### styles.txt ###

Each protocol has an associated plotting style. Those associations are
defined in methods.txt, as described above.

The plotting styles are defined in styles.txt, using three columns:

- style: Numeric ID of style (corresponds to the number given in methods.txt)
- lty: Line type for this style (corresponds to lty parameter in R plot functions)
- pch: Plot symbol for this style (corresponds to pch parameter in R plot functions)

Most scripts will only consider protocols with style ID in the range 1-4. In RGASP,
higher style IDs were used for exploratory purposes and assigned to
protocols not included among those evaluated.

### datasets.txt ###

This files lists the RNA-seq data sets used.

- id: ID for data set. Should consist of word characters only.
- name: Name of dataset. May contain spaces.
- fragments: Total number of read pairs in data set.

Data file locations
-------------------

To be written.

Validation of BAM files
-----------------------

To be continued.