import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns; sns.set()
from scipy.special import factorial 
from scipy.stats import norm, poisson,randint,expon
from sklearn.metrics import classification_report
from Spectrum import Spectrum

# define isotope pattern
def weight(m,k):
    lam = 0.000594 * m - 0.03091
    return np.exp(-lam)*pow(lam,k)/factorial(k)

def iso(x,m,z,sigma):
    result = 0.0
    for k in np.arange(6):
        result += weight(m,k)*norm.pdf(x,loc=(m+k)/z,scale=sigma) # assumes m_neutron = 1       
    return result

def createReferenceBinnedSparse(m,z,NumBinsInit,WindowlengthInit,ResTarget,ResAvg):
    mz = m / z    
    axis = np.linspace(mz-0.5*WindowlengthInit,mz+0.5*WindowlengthInit,num=NumBinsInit)
    intens = iso(axis,m,z,ResAvg)
    
    return Spectrum(axis,intens).bin_to_resolution(ResTarget,min_intensity=1)
        
def createNoiseBinnedSparse(m,z,NumBinsInit,WindowlengthInit,ResTarget,ResAvg):
    mz = m / z    
    axis = np.linspace(mz-0.5*WindowlengthInit,mz+0.5*WindowlengthInit,num=NumBinsInit)
    
    # number peaks
    numPeaks = int(poisson.rvs(mu=1,size=1)[0])+1
    
    # peak locations
    loc = np.sort(randint.rvs(low=0,high=axis.size,size=numPeaks))
    
    # intensities
    intens = expon.rvs(scale=15.,size=numPeaks)
    
    return Spectrum(axis[loc],intens).bin_to_resolution(ResTarget,min_intensity=1)    
    

