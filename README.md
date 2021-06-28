# mzBucket
Locality-sensitive hashing for mass spectrometry data.

## Contents
This repository is structured as follows:
``` bash
.
├── bin
├── data
│   ├── output_create_synthetics
│   ├── output_roc
│   └── output_scalability_study
├── plots
└── src
    ├── cpp
    └── python
        ├── createSynthetics
        ├── meta
        ├── roc
        └── scalability
```
## Build
### Batteries included?
A binary that was compiled on a Ubuntu 20.04 LTS machine with two Intel Xeon Gold 6238 CPUs is provided. 
### Dependencies
The project depends on [eigen](https://gitlab.com/libeigen/eigen) and [opentims](https://github.com/michalsta/opentims).
### Build process
First, clone both projects next to the mzBuild repo.

Then go to the folder opentims++ in the opentims repo and run `make`. This should produce the files zstddeclib.o and sqlite3.o.

Copy both files to /src/cpp/ in the mzBucket repo.

Finally run `make` in  /src/cpp/ in the mzBucket repo.
