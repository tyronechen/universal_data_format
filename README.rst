##########################################################
A flexible, annotation-free biological data representation
##########################################################

Copyright (c) 2020 `Tyrone Chen <https://orcid.org/0000-0002-9207-0385>`_, `Sonika Tyagi <https://orcid.org/0000-0003-0181-6258>`_.

We thank `Yashpal Ramakrishnaiah <https://orcid.org/0000-0002-2213-8348>`_ for suggestions on improving computational efficiency.

Installation
############

With conda::

  git clone "https://gitlab.com/tyagilab/universal_data_format"
  conda create -n <add_your_channel_name_here> \
    pandas tqdm pyfaidx intervaltree \
    bioconductor-biocinstaller r-argparser \
    bioconductor-edger bioconductor-limma bioconductor-rsubread
    bedtools samtools \
    --channel conda-forge --channel bioconda --channel r
  conda activate <add_your_channel_name_here>

Manual install::

  # R packages (within R)
  source("http://bioconductor.org/biocLite.R")
  biocLite(pkgs = c("Rsubread","limma","edgeR"))
  install.packages("argparser")

  # python packages (command line)
  pip install pandas tqdm intervaltree pyfaidx

  # other packages
  # for bedtools, follow the instructions here:
  https://bedtools.readthedocs.io/en/latest/content/installation.html

  # for samtools, follow the instructions here:
  http://www.htslib.org/download/

Usage
#####

Software availability
#####################

All steps, code, parameters, command line arguments and software versions used to generate results in this paper are publicly available in an open source software repository under the `MIT license <https://opensource.org/licenses/MIT>`_, hosted on gitlab at `https://gitlab.com/tyagilab/universal_data_format <https://gitlab.com/tyagilab/universal_data_format>`_. This documentation is provided under a `CC-BY-3.0 AU license <https://creativecommons.org/licenses/by/3.0/au/>`_.

About us
########

`Visit our lab website here <https://bioinformaticslab.erc.monash.edu/>`_. Contact Sonika Tyagi at `sonika.tyagi@monash.edu <mailto:sonika.tyagi@monash.edu>`_.
