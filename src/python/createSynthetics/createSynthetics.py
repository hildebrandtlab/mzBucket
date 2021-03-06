#!/usr/bin/env python3
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns; sns.set()
from tqdm import tqdm
import prototype_lib as lib
from Spectrum import Spectrum
import sys
from itertools import groupby
from functools import reduce

scale = float(sys.argv[1])
numPeaksNoise = int(sys.argv[2])  
outFile = sys.argv[3]

xList = []
yList = []
zList = []
labelList = []
monoList = []
scanList = []

relInt = [0.5,1,0.5]


# labels:
# "noise_1" -- noise in a window without an isotope pattern inside
# "noise_2" -- noise in a window with	 an isotope pattern inside
# "signal"  -- true signal
#

l_n1 = "noise_1"
l_n2 = "noise_2"
l_sig = "signal"

def addLabel(l1,l2):
    return l_sig if ((l1 == l_sig) or (l2 == l_sig)) else l_n2

scan = 0
for m in tqdm(np.arange(150,5000,10)):
    for z in np.arange(1,5):
        if m/z > 150 and m/z < 1999:
            # signal
            #signal = lib.createReferenceBinnedSparse(m,z,1000,10.,1,0.01)
            signal = lib.isoStick(m,z)
          
            # label
            signalLabeled = [ (mz,scale*i,l_sig) for (mz,i) in signal]

            if len(signal) > 0:
                for count in range(3):
                    # modulate signal
                    signalMod = [ (mz,i*relInt[count],l) for (mz,i,l) in signalLabeled]

                    # adding noise 
                    noise = lib.createNoiseBinnedSparse(m,z,1000,10.,1,0.01,numPeaksNoise)
                    
                    # label
                    #noise = [(mz,i,False) for (mz,i) in noise]
                    noise = [(mz,i,l_n2) for (mz,i) in noise]
                    # concat lists
                    signalMod += noise
                    
                    # group by 
                    grouped = []
                    
                    for mzKey, values in groupby(signalMod, lambda x: int(1000*x[0])):
                        listVal = list(values)                                             
                        grouped.append(\
                                (float(mzKey/1000),\
                                reduce(lambda i1, i2: i1 + i2, [ i for (mz,i,l) in listVal]),\
                                reduce(lambda lab1, lab2: addLabel(lab1,lab2), [ l for (mz,i,l) in listVal]))\
                                )
                    
                    # sort by
                    sortedTup = sorted(grouped, key=lambda tup: tup[0])
                    
                    # append on lists
                    for (mz,i,label) in sortedTup:
                        if mz > 150 and mz < 1999:
                           xList.append(mz)
                           yList.append(i)
                           zList.append(z)
                           monoList.append(m/z)
                           labelList.append(label)
                           scanList.append(scan)

                    scan +=1 
            for count in range(10):
                # noise 
                noise = lib.createNoiseBinnedSparse(m,z,1000,10.,1,0.01,numPeaksNoise)
                for (mz,i) in noise:
                    if mz > 150 and mz < 1999:
                      xList.append(mz)
                      yList.append(i)
                      zList.append(z)
                      monoList.append(m/z)
                      labelList.append(l_n1)
                      scanList.append(scan)

                scan +=1 


df_reference = pd.DataFrame.from_dict({'scan':scanList,'mz':xList,'i':yList,'label':labelList})

df_reference.to_csv(outFile,index=False)

