#!/usr/bin/env python3

import subprocess as sp

scale = 2000
labelNoise = "False"
numPeaksNoise = 4 

for scale in [1000,500,250,125,64,32]:
	print("Running scale {}".format(scale))
	
	syntheticsFile =  "../../../data/synthetics_{}_{}_{}.csv".format(scale,labelNoise,numPeaksNoise)

	plotFile =  "../../../plots/ROC_PP_{}_{}_{}.pdf".format(scale,labelNoise,numPeaksNoise)
	# create synths
	sp.run(["../createSynthetics/createSynthetics.py",str(scale),str(labelNoise),str(numPeaksNoise),syntheticsFile])

	# do hashing
	sp.run(["../roc/SNRonSynthetics.py",syntheticsFile,plotFile])
