# Imports ---------------------------------------
import sys
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'  # or any {'0', '1', '2'}

import tensorflow as tf
import numpy as np

from proteolizarddata.data import PyTimsDataHandle, MzSpectrum, TimsFrame
from proteolizardalgo.hashing import TimsHasher

import pandas as pd
from sklearn.metrics import confusion_matrix, roc_auc_score,roc_curve
from scipy.stats import median_abs_deviation

from tqdm import tqdm
from matplotlib import pyplot as plt
import seaborn as sns
from adjustText import adjust_text

# Plot settings ---------------------------------------
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


# Parse arguments ---------------------------------------
inFile = sys.argv[1]
outFile = sys.argv[2]

# Load synthetics ---------------------------------------
df = pd.read_csv(inFile)
df['label_bin'] = df.apply(lambda r: r['label'] != 'noise_1',axis=1)

# Binning 
window_length = 10.
df['bins'] = np.floor(df.mz/window_length).astype(np.int32)
df['id'] = df.index.values

# Rounding
syntheticsFrame = TimsFrame(None,-1, df.scan.values, df.mz.values, df.i.values.astype(np.int32), np.zeros_like(df.i,dtype=np.int32), np.zeros_like(df.i))
s, b, F = syntheticsFrame.get_dense_windows(min_intensity=0, resolution=1, min_peaks=1, overlapping=False)

# Hashing ---------------------------------------

def get_roc_point(num_hashes,len_single_hash,s,b,F,df):

    # generate hasher and get random matrix
    hasher = TimsHasher(num_hashes, len_single_hash, 13, 1, num_dalton=10)

    # get hash keys
    H = hasher.calculate_keys(F)

    # get collisions
    ret_bins, ret_scans = hasher.hash_ptr.calculateCollisions(H, s, b)
    df_res = pd.DataFrame({'ret_bins':ret_bins,'ret_scans':ret_scans})


    ## Join results
    merged = pd.merge(left=df,right=df_res,right_on=["ret_scans","ret_bins"],left_on=["scan","bins"],indicator=True)
    called_peak_ids = set(merged.id.values)

    df['label_hat'] = df.apply(lambda r: r['id'] in called_peak_ids,axis=1)
    tn, fp, fn, tp = confusion_matrix(df.label_bin.values, df.label_hat.values).ravel()
    tpr = tp / (tp + fn)
    fpr = fp / (fp + tn)
    
    return fpr,tpr

ors_list = [64,32,16]
ands_list = [16,32,64]

# store results
kList = []
lList = []
tprList = []
fprList = []

for ors in tqdm(ors_list):
    for ands in ands_list:
        fpr,tpr = get_roc_point(ors,ands,s,b,F,df)

        kList.append(ors)
        lList.append(ands)
        tprList.append(tpr)
        fprList.append(fpr)

dfLSH = pd.DataFrame.from_dict({'ands':lList,'ors':kList,'tpr':tprList,'fpr':fprList})

# Thresholding ---------------------------------------

def getNoiseEstimate(yList):
    if len(yList) > 1:
        sigma = median_abs_deviation(yList)
    else:
        sigma = np.inf

    return sigma

def get_roc_thres(df_arg):
    # copy
    df = df_arg.copy()
    
    # Calculate noise estimate on windows
    windowDF = df[['scan','i']].groupby('scan',as_index=False).agg(list)
    windowDF['sigma'] = windowDF.apply(lambda r: getNoiseEstimate(r['i']),axis=1)
    
    # store sigma for each window
    dictSigma = {k:v for (k,v) in zip(list(windowDF.scan.values),list(windowDF.sigma.values))}
    
    # add sigma and SNR for each row in df, i.e. per peak
    df['sigma'] = df.apply(lambda r:dictSigma[r['scan']],axis=1) 
    df['SNR'] = df['i'] / df['sigma']
   
    # remove pathological values    
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    df = df.dropna(subset=["SNR"], how="all")
    
    # set labels
    df['yTrue'] = df.apply(lambda r:1 if r['label'] == "signal" else 0,axis=1)
    
    # do not consider noise peaks within a pattern
    df = df[df['label'] != "noise_2"]
      
    fpr, tpr, thresholds = roc_curve(y_true=df.yTrue.values,y_score=df.SNR.values,pos_label=1)
    return (fpr, tpr, thresholds)


res = get_roc_thres(df)

# Plot results ---------------------------------------
plt.xlim([0,1])
plt.ylim([0,1])
plt.plot(np.linspace(0,1),np.linspace(0,1),'--',alpha=0.5)
plt.xlabel("False postive rate",size=MEDIUM_SIZE)
plt.ylabel("True positive rate",size=MEDIUM_SIZE)
plt.title("Classification results for different $m$ and $n$.",size=BIGGER_SIZE);

# LSH
tprVals = dfLSH['tpr'].values
fprVals = dfLSH['fpr'].values
keys = [(l,k) for (l,k) in zip( dfLSH['ands'].values, dfLSH['ors'].values)]
plt.scatter(fprVals,tprVals,s=5,label="LSH")

# Text 
texts = [plt.text(fprVals[i],tprVals[i],k) for i,k in enumerate(keys)]
adjust_text(texts, arrowprops=dict(arrowstyle='->', color='red'))

# SNR
thres_fpr, thres_tpr,_ = res
plt.plot(thres_fpr,thres_tpr,label="SNR",color="C1")

plt.legend(loc='lower right')
plt.tight_layout()
        
plt.savefig(outFile,dpi=dpi)
