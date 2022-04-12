# mzBucket
Locality-sensitive hashing for mass spectrometry data.

## Contents
This repository is structured as follows:
``` bash
.
├── data
│   ├── output_create_synthetics
│   ├── output_roc
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


### Dependencies
The project depends on [eigen](https://gitlab.com/libeigen/eigen) and [opentims](https://github.com/michalsta/opentims).
### Build process
First, clone both projects next to the mzBuild repo.

Then go to the folder opentims++ in the opentims repo and run `make`. This should produce the files zstddeclib.o and sqlite3.o.

Copy both files to /src/cpp/ in the mzBucket repo.

Finally run `make` in  /src/cpp/ in the mzBucket repo.

## Data 
The data set used is available for download [here](https://zenodo.org/record/5036526#.YlWE91xBw5l). Please unzip it into the data folder to ensure reproducible runs of the scripts.
