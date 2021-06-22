#!/usr/bin/env python3

import subprocess as sp

numPeaksNoise = 4 

for scale in [1000,500,250,125,64,32]:
	print("Running scale {}".format(scale))
	syntheticsFile =  "../../../data/output_create_synthetics/synthetics_{}_{}.csv".format(scale,numPeaksNoise)
	rocFile =  "../../../data/output_roc/roc_{}_{}.csv".format(scale,numPeaksNoise)

	# create synths
	#sp.run(["../createSynthetics/createSynthetics.py",str(scale),str(numPeaksNoise),syntheticsFile])

	# do hashing
	sp.run(["../roc/LSHonSynthetics.py",syntheticsFile,rocFile])

	# do plotting 
	sp.run(["../roc/plotROC.py",rocFile,"../../../plots/ROC_{}_{}.pdf".format(scale,numPeaksNoise)])

