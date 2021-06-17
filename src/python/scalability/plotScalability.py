from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import sys

inFile = sys.argv[1]
outFile = sys.argv[2]

df = pd.read_csv(inFile)
#df['k'] = 5
#df['l'] = 6
df = df.assign(Amplification=df.apply(lambda r: "({},{})".format(int(r['l']),int(r['k'])),axis=1))
df["Amplification"] = df["Amplification"].astype("category")
df = df.sort_values('k')
print(df)

widthMM = 170
widthInch = widthMM / 25.4
ratio = 0.66666
heigthInch = ratio*widthInch

SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 12

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
sns.set_style("ticks")

dpi = 300


fig, ax = plt.subplots(figsize=(widthInch, heigthInch), dpi=dpi, facecolor='w', edgecolor='k')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_yscale('log')
ax.set_xscale('log')
plt.title("Scalability for different $m$ and $n$.",size=BIGGER_SIZE);

sns.scatterplot(data=df,x='numThreads',y='runTimeSec',hue='Amplification',style='Amplification')#,hue_order=["(30,22)","(30,32)","(30,64)"])
#sns.scatterplot(data=df,x='numThreads',y='runTime',hue='hue',style='hue')
plt.xlabel("Number of threads used",size=MEDIUM_SIZE)
plt.ylabel("Runtime [s]",size=MEDIUM_SIZE)

plt.legend()
plt.tight_layout()
plt.savefig(outFile,dpi=dpi)

