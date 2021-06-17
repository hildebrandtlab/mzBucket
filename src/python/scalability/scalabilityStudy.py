import subprocess
import time
import pandas as pd
import numpy as np
import sys

def run(numThreads,k,l):
    start = time.time()
    # ./mzLSH -f 1009 -k 30 -l 30 -d /share/massSpec/bruker/dda/M210115_001_Slot1-1_1_850.d/ -b ~/timsj/src/main/resources/libtimsdata.so -t 32
    runCommand = ["../mzLSH","-f","1009","-k","{}".format(k),"-l","{}".format(l),"-d","/share/massSpec/bruker/dda/M210115_001_Slot1-1_1_850.d/","-b","/home/bob/timsj/src/main/resources/libtimsdata.so","-t"]
    runCommand.append("{}".format(numThreads))
    print(runCommand)
    subprocess.run(runCommand)
    end = time.time()
    return end - start

# outFile
outfile = sys.argv[1]

with open("{}_tmp".format(outfile), 'a') as the_file:
    the_file.write('numThreads,runTime,l,k\n')

# run parameters
numThreads = [2,4,8,16,32,64,88]
ampParameters = [(30,64),(30,30),(30,22)]

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


