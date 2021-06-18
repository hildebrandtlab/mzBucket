//
// Created by David Teschner on 07.06.21.
//

#ifndef TIMSLSH_HASHING_H
#define TIMSLSH_HASHING_H

#include <iostream>
#include <cassert>
#include <random>
//#include "/home/administrator/Documents/promotion/eigen/Eigen/Dense"
//#include "/home/administrator/Documents/promotion/eigen/Eigen/Sparse"

//#include "/Users/davidteschner/Documents/promotion/eigen/Eigen/Dense"
//#include "/Users/davidteschner/Documents/promotion/eigen/Eigen/Sparse"

#include "/home/dateschn/eigen/Eigen/Dense"
#include "/home/dateschn/eigen/Eigen/Sparse"

#define assertm(exp, msg) assert(((void)msg, exp))

/**
 * create eigen sparse vector from vectorized mz spectrum
 * @param mzVector : vectorized mz spectrum to convert
 * @param numRows : dimensionality of vector
 * @return : a sparse eigen vector suited for fast vectorized operations
 */
Eigen::SparseVector<double> toSparseVector(const MzVector& mzVector, int numRows){

    auto sparseVec = Eigen::SparseMatrix<double>(numRows, 1);

    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(mzVector.indices.size());

    for(int i = 0; i < mzVector.indices.size(); i++)
        tripletList.emplace_back(mzVector.indices[i], 0, mzVector.values[i]);

    sparseVec.setFromTriplets(tripletList.begin(), tripletList.end());

    return sparseVec;
}

/**
 * create a string representation of bits for easy hashing
 * @param boolVector vector to create bit string from
 * @return bit string
 */
std::string boolVectorToString(const std::vector<bool> &boolVector, int bin, bool restricted){
    std::string ret = "";
    for(const auto& b: boolVector)
        b ? ret.append("1") : ret.append("0");

    // this is a hard restriction to its on mass bin only for collision
    if(restricted)
        ret.append(std::to_string(bin));
    // this is a soft restriction to all windows with same offset
    else
        bin > 0 ? ret.append("1") : ret.append("0");
    return ret;
}


/**
 * calculate integer keys from vectors of boolean
 * @param hashes set of bool vectors representing the result from lsh probing
 * @return a set of ints representing the keyspace of a given mz spectrum or window
 */
std::vector<int> calculateKeys(const std::vector<std::vector<bool>> &hashes, int bin, bool restricted){
    std::vector<int> retVec;
    retVec.reserve(hashes.size());
    for(const auto& h: hashes){
        // const auto bitString = boolVectorToString(h, bin, restricted);
        int hash = std::hash<std::string>{}(boolVectorToString(h, bin, restricted));
        retVec.push_back(hash);
    }
    return retVec;
}

/**
 *
 * @param sparseSpectrumVector
 * @param M random matrix of shape ((k * l) * mzSpace)
 * @param k number of ANDs
 * @param l number of ORs
 * @return k vectors of l bits
 */
std::vector<std::vector<bool>> calculateSignumVector(const Eigen::SparseVector<double>& sparseSpectrumVector, const Eigen::MatrixXd& M,
                                                     int k, int l){

    // check for compatible settings
    assertm(k * l == M.rows(), "dimensions of random vector matrix and banding differ.");

    std::vector<std::vector<bool>> retVec;
    retVec.reserve(k);

    // heavy lifting happens here, calculate dot products between random vectors and spectrum
    auto r = M * sparseSpectrumVector;

    // calculate signs from dot products
    std::vector<bool> bVec;
    bVec.reserve(l);

    for(int i = 0; i < r.size(); i++)
        bVec.push_back(r[i] > 0);

    // rest of the code is for grouping of results only
    std::vector<bool> leVec;
    leVec.reserve(l);

    for(int i = 0; i < bVec.size(); i++){

        if(i % l == 0 && i != 0){
            retVec.push_back(leVec);
            leVec.clear();
        }
        leVec.push_back(bVec[i]);
    }
    retVec.push_back(leVec);

    return retVec;
}

/**
 * check if key is already in map
 * @param hashKeys map to check
 * @param key to look for
 * @return true if key is in map else false
 */
bool hasValue(const std::map<int, bool>& hashKeys, const int key){
    return hashKeys.count(key);
}

/**
 * check if key is already in map
 * @param containsBins map to check
 * @param key to look for
 * @return true if key is in map else false
 */
bool hasValue(const std::map<int, std::set<int>>& containsBins, const int key){
    return containsBins.count(key);
}

/**
 * check if a given key occurred in map more then once
 * @param keys
 * @param containsMap
 * @return
 */
bool containsKeyMultiple(const std::vector<int>& keys, const std::map<int, bool>& containsMap){
    bool moreThenOnce = false;
    for(const auto key: keys){
        if(containsMap.at(key)){
            moreThenOnce = true;
        }
    }
    return moreThenOnce;
}


/**
 *
 * @param specMap
 * @param numThreads
 * @param windowLength
 * @param overlapping
 * @param k
 * @param l
 */
std::map<int, std::set<int>> getHashes(const std::map<int, MzSpectrum>& specMap, int numThreads, double windowLength,
                                       bool overlapping, int k, int l, bool normalize, bool sqr, bool restricted, bool verbose) {
    std::vector<std::pair<int, MzSpectrum>> mzSpectra;
    mzSpectra.reserve(specMap.size());

    std::vector<std::pair<int, MzSpectrum>> mzSpectraBinned;
    mzSpectraBinned.reserve(specMap.size());
    mzSpectraBinned.resize(specMap.size());

    std::vector<std::pair<int, std::map<int, MzSpectrum>>> binnedWindows;
    binnedWindows.reserve(specMap.size());
    binnedWindows.resize(specMap.size());

    std::vector<std::pair<int, std::map<int, MzVector>>> vectorizedWindows;
    vectorizedWindows.reserve(specMap.size());
    vectorizedWindows.resize(specMap.size());

    std::vector<std::pair<int, std::map<int, Eigen::SparseVector<double>>>> vectorizedWindowsEigen;
    vectorizedWindowsEigen.reserve(specMap.size());
    vectorizedWindowsEigen.resize(specMap.size());

    std::vector<std::pair<int, std::map<int, std::vector<int>>>> hashKeys;
    hashKeys.reserve(specMap.size());
    hashKeys.resize(specMap.size());

    // create random number generator
    std::default_random_engine generator(42);
    std::normal_distribution<double> distribution(0, 1);
    auto normal = [&] (int) {return distribution(generator);};

    // fill matrix with random numbers
    const Eigen::MatrixXd M = Eigen::MatrixXd::NullaryExpr(k * l, 20000, normal);

    for(const auto& [key, value]: specMap){
        mzSpectra.push_back({key, value});
    }

    if(verbose)
        std::cout << "Number of spectra: " << mzSpectra.size() << std::endl;

    // binning
    #pragma omp parallel for num_threads(numThreads)
    for(int i = 0; i < specMap.size(); i++){
        MzSpectrum binnedSpec = binToResolution(mzSpectra[i].second, 1);
        mzSpectraBinned[i] = {mzSpectra[i].first, binnedSpec};
    }

    // windows generation
    #pragma omp parallel for num_threads(numThreads)
    for(int i = 0; i < specMap.size(); i++){
        auto windowCollection = groupToWindows(mzSpectraBinned[i].second, windowLength, overlapping);
        binnedWindows[i] = {mzSpectraBinned[i].first, windowCollection};
    }

    int windowCounter = 0;
    for(const auto& [scan, windows]: binnedWindows){
        windowCounter += windows.size();
    }

    if(verbose)
        std::cout << "Number of windows: " << windowCounter << std::endl;

    // window vectorization
    #pragma omp parallel for num_threads(numThreads)
    for(int i = 0; i < specMap.size(); i++){

        std::map<int, MzVector> tmpVec;
        for(const auto& [first, second]: binnedWindows[i].second){
            tmpVec[first] = toVectorizedBinnedSpectrum(second, 1, normalize, sqr);
        }

        vectorizedWindows[i] = {binnedWindows[i].first, tmpVec};
    }

    // window to eigen vec
    #pragma omp parallel for num_threads(numThreads)
    for(int i = 0; i < specMap.size(); i++){

        std::map<int, Eigen::SparseVector<double>> tmpVec;
        for(const auto& [first, second]: vectorizedWindows[i].second){
            tmpVec[first] = toSparseVector(second, 20000);
        }

        vectorizedWindowsEigen[i] = {vectorizedWindows[i].first, tmpVec};
    }

    if(verbose){
        std::cout << "______________________________" << std::endl;
        std::cout << "" << std::endl;

        std::cout << "Starting LSH." << std::endl;
        std::cout << "______________________________" << std::endl;
    }

    // key generation
    #pragma omp parallel for num_threads(numThreads)
    for(int i = 0; i < specMap.size(); i++){

        std::map<int, std::vector<int>> tmpVec;
        for(const auto& [first, second]: vectorizedWindowsEigen[i].second){
            tmpVec[first] = calculateKeys(calculateSignumVector(second, M, k, l), first, restricted);
        }

        hashKeys[i] = {vectorizedWindowsEigen[i].first, tmpVec};
    }

    int keyCounter = 0;
    std::set<int> keySet;

    for(const auto& [scan, windows]: hashKeys){
        for(const auto& [bin, keys]: windows){
            keyCounter += keys.size();
            for(const auto& key: keys)
                keySet.insert(key);
        }
    }

    if(verbose){
        std::cout << "Number of keys       : " << keyCounter << std::endl;
        std::cout << "Number of unique keys: " << keySet.size() << std::endl;
    }

    std::map<int, bool> filterMap;

    for(const auto& [key, values]: hashKeys){
        for(const auto& [k, v]: values){
            for(const auto& hashKey: v){
                if(hasValue(filterMap, hashKey)){
                    filterMap[hashKey] = true;
                } else {
                    filterMap[hashKey] = false;
                }
            }
        }
    }

    std::map<int, std::set<int>> checkSet;

    for(const auto&[scan, windows]: hashKeys){
        for(const auto&[bin, keys]: windows){
            if(containsKeyMultiple(keys, filterMap)){
                auto s = checkSet[scan].insert(bin);
            }
        }
    }

    return checkSet;
}
#endif //TIMSLSH_HASHING_H
