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

outFile = sys.argv[1]

xList = []
yList = []
zList = []
labelList = []
monoList = []
scanList = []

relInt = [0.5,1,0.5]

scan = 0
for m in tqdm(np.arange(150,5000,10)):
    for z in np.arange(1,5):
        if m/z > 150 and m/z < 1999:
            # signal
            tuples = lib.createReferenceBinnedSparse(m,z,1000,10.,1,0.01)
            
            # label
            tuples = [ (mz,i,True) for (mz,i) in tuples]

            if len(tuples) > 0:
                for count in range(3):
                    # modulate signal
                    tuples = [ (mz,i*relInt[count],l) for (mz,i,l) in tuples]

                    # adding noise 
                    noise = lib.createNoiseBinnedSparse(m,z,1000,10.,1,0.01)
                    
                    # label
                    noise = [(mz,i,False) for (mz,i) in noise]
                    #noise = [(mz,i,True) for (mz,i) in noise]
                    # concat lists
                    tuples += noise
                    
                    # group by 
                    grouped = []
                    
                    for mzKey, values in groupby(tuples, lambda x: x[0]):
                        listVal = list(values)                                             
                        grouped.append(\
                                (mzKey,\
                                reduce(lambda i1, i2: i1 + i2, [ i for (mz,i,l) in listVal]),\
                                reduce(lambda lab1, lab2: lab1 or lab2, [ l for (mz,i,l) in listVal]))\
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
                tuples = lib.createNoiseBinnedSparse(m,z,1000,10.,1,0.01)
                for (mz,i) in tuples:
                    if mz > 150 and mz < 1999:
                      xList.append(mz)
                      yList.append(i)
                      zList.append(z)
                      monoList.append(m/z)
                      labelList.append(False)
                      scanList.append(scan)

                scan +=1 


df_reference = pd.DataFrame.from_dict({'scan':scanList,'mz':xList,'i':yList,'label':labelList})

df_reference.to_csv(outFile,index=False)

