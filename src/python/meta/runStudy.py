#!/usr/bin/env python3

import subprocess as sp

scale = 2000
labelNoise = "True"
numPeaksNoise = 1


syntheticsFile =  "../../../data/synthetics_{}_{}_{}.csv".format(scale,labelNoise,numPeaksNoise)
rocFile =  "../../../data/roc_{}_{}_{}.csv".format(scale,labelNoise,numPeaksNoise)

# create synths
sp.run(["../createSynthetics/createSynthetics.py",str(scale),str(labelNoise),str(numPeaksNoise),syntheticsFile])

# do hashing
sp.run(["../roc/LSHonSynthetics.py",syntheticsFile,rocFile])

# do plotting 
sp.run(["../roc/plotROC.py",rocFile,"../../../plots/ROC_{}_{}_{}.pdf".format(scale,labelNoise,numPeaksNoise)])

