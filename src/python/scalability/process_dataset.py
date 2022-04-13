#!/usr/bin/env python3

import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'  # or any {'0', '1', '2'}
# disable gpu
os.environ["CUDA_VISIBLE_DEVICES"] = "-1"

import sys
from tqdm import tqdm
import numpy as np

from proteolizarddata.data import PyTimsDataHandle, MzSpectrum
from proteolizardalgo.hashing import TimsHasher

# Parse args --------------------------
if len(sys.argv) != 4:
    print("Usage python process_dataset.py inFile numHashes len_single_hash")
    sys.exit(1)
    
inFile = sys.argv[1]

# set hashing hyper parameters
num_hashes = int(sys.argv[2])
len_single_hash = int(sys.argv[3])

# Setup  --------------------------
# create data handle and load a precursor frame 
dh = PyTimsDataHandle(inFile)

# generate hasher and get random matrix
hasher = TimsHasher(num_hashes, len_single_hash, 13, 1, num_dalton=10)

frame_list = [dh.get_frame(i) for i in np.random.choice(dh.precursor_frames,100)] # study four has 400 frames

# Run  --------------------------
for frame in tqdm(frame_list):
    # get data tensors
    s, b, F = frame.get_dense_windows(min_intensity=30, resolution=1, min_peaks=3, overlapping=False)
    
    # get hash keys
    H = hasher.calculate_keys(F)
    
    # get collisions
    ret_bins, ret_scans = hasher.calculate_collisions(H, s, b)