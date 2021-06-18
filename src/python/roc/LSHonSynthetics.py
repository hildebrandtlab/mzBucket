import subprocess
import time
import pandas as pd
import numpy as np
import sys

def run(l,k):
    start = time.time()
    # ./mzBucket -f 1 -k 20 -l 20 -c synthetics.csv -t 32
    runCommand = ["/share/massSpec/syntheticsLSH/mzBucket","-f","1","-k","{}".format(k),"-l","{}".format(l),"-c","/share/massSpec/syntheticsLSH/synthetics.csv","-t","16","-v","false"]
    print("// {}".format(runCommand))
    res = subprocess.run(runCommand)
    return res

# outFile
outfile = sys.argv[1]

# run parameters
ampParameters = [(30,22)]

# store results
runTime = []
kList = []
lList = []
threadsUsed = []
for (l,k) in ampParameters:
    duration = run(l,k) 
        #print("l {} k {} time {}".format(l,k,duration))
        #runTime.append(duration)
        #kList.append(k)
        #lList.append(l)
        #threadsUsed.append(nt)
        #with open("{}_tmp".format(outfile), 'a') as the_file:
        #    the_file.write('{},{},{},{}\n'.format(nt,duration,l,k))

#df = pd.DataFrame.from_dict({'numThreads':threadsUsed,'runTimeSec':runTime,'k':kList,'l':lList})
#df.to_csv(outfile)
#print(df)


