import pandas as pd
import numpy as np
import sys
from matplotlib import pyplot as plt
from scipy.stats import median_abs_deviation
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve
import seaborn as sns

def getNoiseEstimate(yList):
    if len(yList) > 1:
        sigma = median_abs_deviation(yList)
    else:
        sigma = np.inf

    return sigma

inFile = sys.argv[1]
outFile = sys.argv[2]

df = pd.read_csv(inFile)
print(df)

# Calculate noise estimate on windows
windowDF = df[['scan','intensity']].groupby('scan',as_index=False).agg(list)
print(windowDF)

windowDF['sigma'] = windowDF.apply(lambda r: getNoiseEstimate(r['intensity']),axis=1)
print(windowDF)

dictSigma = {k:v for (k,v) in zip(list(windowDF.scan.values),list(windowDF.sigma.values))}

df['sigma'] = df.apply(lambda r:dictSigma[r['scan']],axis=1)
df['yTrue'] = df.apply(lambda r:1 if r['label'] else 0,axis=1)
df['SNR'] = df['intensity'] / df['sigma']

print(df)

fpr, tpr, thresholds = roc_curve(y_true=df.yTrue.values,y_score=df.SNR.values)

dfROC = pd.DataFrame.from_dict({'tpr':tpr,'fpr':fpr,'t':thresholds})
print(dfROC)

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
#ax.set_yscale('log')
#ax.set_xscale('log')
#plt.title("Scalability for different $m$ and $n$.",size=BIGGER_SIZE);
plt.step(tpr,fpr)

plt.xlim([0,1])
plt.ylim([0,1])
plt.plot(np.linspace(0,1),np.linspace(0,1),'--',alpha=0.5)
plt.xlabel("False postive rate",size=MEDIUM_SIZE)
plt.ylabel("True positive rate",size=MEDIUM_SIZE)

#plt.legend()
plt.tight_layout()
plt.savefig(outFile,dpi=dpi)

