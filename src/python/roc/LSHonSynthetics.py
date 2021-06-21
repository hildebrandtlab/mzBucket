#!/usr/bin/env python3
import subprocess as sp
import time
import pandas as pd
import numpy as np
import sys
import re
from tqdm import tqdm

def run(l,k,inFile):
    start = time.time()
    # ./mzBucket -f 1 -k 20 -l 20 -c synthetics.csv -t 32
    runCommand = ["../../../bin/mzBucket","-f","1","-k","{}".format(k),"-l","{}".format(l),"-c","{}".format(inFile),"-t","32","-v","false","-r","false","-w","10"]
    #print("// {}".format(runCommand))
    res = str(sp.check_output(runCommand))
    #print("res {}".format(res))
    counts = re.findall(r'\d+', res)
    #print("counts {}".format(counts))
    tp = int(counts[0])
    tn = int(counts[1])
    fp = int(counts[2])
    fn = int(counts[3])
    return (tp,fp,tn,fn)

# inFile 
inFile = sys.argv[1]

# outFile
outfile = sys.argv[2]

# run parameters
ampParameters = [(30, 25),(30, 28),(30, 30),(30, 22),(30, 32),(30, 38),(30, 45),(30, 64)]
#ampParameters = [(30,25)]

# store results
kList = []
lList = []
tprList = []
fprList = []
for (l,k) in tqdm(ampParameters):
    (tp,fp,tn,fn) = run(l,k,inFile) 
 
    tpr = tp / (tp + fn)
    fpr = fp / (fp + tn)

    kList.append(k)
    lList.append(l)
    tprList.append(tpr)
    fprList.append(fpr)

    with open("{}_tmp".format(outfile), 'a') as the_file:
            the_file.write('{},{},{},{}\n'.format(l,k,tpr,fpr))

df = pd.DataFrame.from_dict({'l':lList,'k':kList,'tpr':tprList,'fpr':fprList})
df.to_csv(outfile,index=False)
print(df)


