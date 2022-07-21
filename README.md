# mzBucket
Locality-sensitive hashing for mass spectrometry data.

## Publication
See [here](https://doi.org/10.1186/s12859-022-04833-5) for the corresponding publication.

## Contents
This repository is structured as follows:
``` bash
.
├── data
│   ├── output_create_synthetics
│   └── output_scalability_study
├── plots
└── src
    └── python
        ├── createSynthetics
        ├── meta
        ├── roc
        └── scalability
```
## Build
The implementation has been updated and shifted to the proteolizard framework, more precisely the components [pyproteolizard-data](https://github.com/theGreatHerrLebert/proteolizard-data) for data-access and [pyproteolizard-algorithm](https://github.com/theGreatHerrLebert/proteolizard-algorithm) for the LSH. See in the respective repositories for buils instructions. 

## Data 
The data set used is available for download [here](https://zenodo.org/record/5036526#.YlWE91xBw5l). Please unzip it into the data folder to ensure reproducible runs of the scripts.
