import subprocess
import time
import pandas as pd
import numpy as np
import sys

def run(numThreads,l,k):
    percent = numThreads*100
    start = time.time()
    #cpulimit -f -l 1000 -- python process_dataset.py ../../data/M210115_001_Slot1-1_1_850.d/ l k 
    runCommand = ["cpulimit","-f","-l","{}".format(percent),"--","python","process_dataset.py","../../data/M210115_001_Slot1-1_1_850.d/","{}".format(l),"{}".format(k)]

    print(runCommand)
    subprocess.run(runCommand)
    end = time.time()
    return end - start

# outFile
outfile = sys.argv[1]

with open("{}_tmp".format(outfile), 'a') as the_file:
    the_file.write('numThreads,runTime,l,k\n')

# run parameters
numThreads = [2,4,8,16,32,64,128]
#ampParameters = [(30,64),(30,30),(30,22)] -> up to study #4
ampParameters = [(32,64),(32,32),(16,32)]

# store results
runTime = []
kList = []
lList = []
threadsUsed = []
for nt in np.flip(numThreads):
    for (l,k) in ampParameters:
        duration = run(nt,l,k) 
        print("l {} k {} time {}".format(l,k,duration))
        runTime.append(duration)
        kList.append(k)
        lList.append(l)
        threadsUsed.append(nt)
        with open("{}_tmp".format(outfile), 'a') as the_file:
            the_file.write('{},{},{},{}\n'.format(nt,duration,l,k))

df = pd.DataFrame.from_dict({'numThreads':threadsUsed,'runTimeSec':runTime,'k':kList,'l':lList})
df.to_csv(outfile)
print(df)