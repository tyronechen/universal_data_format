# A flexible, annotation-free biological data representation

> NOTE: [The main repository is on gitlab](https://gitlab.com/tyagilab/universal_data_format). It is also [mirrored on github](https://github.com/tyronechen/universal_data_format) but please submit any issues to the main gitlab repository only.

Copyright (c) 2020-2021 <a href="https://orcid.org/0000-0002-9207-0385">Tyrone Chen <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>, <a href="https://orcid.org/0000-0003-0181-6258">Sonika Tyagi <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>

Chen T and Tyagi S. An annotation-free format for representing multimodal data features  [version 1; not peer reviewed]. F1000Research 2021, 10(ISCB Comm J):657 (poster) [(https://doi.org/10.7490/f1000research.1118642.1)](https://doi.org/10.7490/f1000research.1118642.1) 

## Contents

[[_TOC_]]

## Acknowledgements

We thank <a href="https://orcid.org/0000-0002-2213-8348">Yashpal Ramakrishnaiah <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a> for suggestions on improving computational efficiency, and <a href="https://orcid.org/0000-0002-1509-4981">Matúš Kalaš <img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a> for valuable suggestions and feedback.

## Software availability

All steps, code, parameters, command line arguments and software versions used to generate results will be made available upon publication in an open source software repository under the [MIT license](https://opensource.org/licenses/MIT), hosted on gitlab at [https://gitlab.com/tyagilab/universal_data_format](https://gitlab.com/tyagilab/universal_data_format). This documentation is provided under a [CC-BY-3.0 AU license](https://creativecommons.org/licenses/by/3.0/au/).

## About us

[Visit our lab website here](https://bioinformaticslab.erc.monash.edu/). Contact Sonika Tyagi at [sonika.tyagi@monash.edu](mailto:sonika.tyagi@monash.edu).

## Installation

### Quick install

With docker and singularity:

*To be written*

### Conda install

With conda:

```
git clone "https://gitlab.com/tyagilab/universal_data_format"
conda create -n <add_your_channel_name_here> \
  pandas tqdm pyfaidx intervaltree \
  bioconductor-biocinstaller r-argparser \
  bioconductor-edger bioconductor-limma bioconductor-rsubread
  bedtools samtools \
  --channel conda-forge --channel bioconda --channel r
conda activate <add_your_channel_name_here>
```

### Manual install

Manual install:

```
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
```

## Usage

*To be written*
