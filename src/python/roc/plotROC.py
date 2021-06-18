import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns; sns.set()
from sklearn.metrics import roc_curve, auc, RocCurveDisplay, roc_auc_score
import sys

# Plot settings --------------------------------
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

# -----------------------------------------------

inFile = sys.argv[1]
outFile = sys.argv[2]

dfLSH = pd.read_csv(inFile)
print(dfLSH)

tprVals = dfLSH['tpr'].values
fprVals = dfLSH['fpr'].values
keys = [(l,k) for (l,k) in zip( dfLSH['l'].values, dfLSH['k'].values)]

fig, ax = plt.subplots(figsize=(widthInch, heigthInch), dpi=dpi, facecolor='w', edgecolor='k')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.scatter(fprVals,tprVals)
plt.xlim([0,1])
plt.ylim([0,1])
plt.plot(np.linspace(0,1),np.linspace(0,1),'--',alpha=0.5)
plt.xlabel("False postive rate",size=MEDIUM_SIZE)
plt.ylabel("True positive rate",size=MEDIUM_SIZE)
plt.title("Classification results for different $m$ and $n$.",size=BIGGER_SIZE);
#plt.text(0.8,0.05,"AUC = {0:.3g}".format(auc(*zip(*sorted(list(zip(fprVals,tprVals)),key=lambda t: t[0])))),size=BIGGER_SIZE)
for i,k in enumerate(keys):
        plt.text(fprVals[i],tprVals[i],k)
        plt.tight_layout()
        plt.savefig(outFile,dpi=dpi)

