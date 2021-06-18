//
// Created by David Teschner on 07.06.21.
//

#ifndef TIMSLSH_SPECTRUM_H
#define TIMSLSH_SPECTRUM_H

#include <math.h>

/**
 * container for timsTOF mz spectra
 */
struct MzSpectrum {
    int frameId {}; // rt coordinate
    int scanId {}; // dt coordinate
    std::vector<double> mz; // vector of mz values
    std::vector<double> intensity; // vector of associated intensities
};

/**
 * container for vectorized timsTOF mz spectra
 */
struct MzVector {
    int resolution; // number of decimals spectrum was binned to
    int frameId; // rt coordinate
    int scanId; // dt coordinate
    std::vector<int> indices; // mz indices, can be seen as vector entries
    std::vector<double> values; // vector of associated intensities
};

/**
 * check if a given windowId is already present in map
 * @param windowCollection collection to check
 * @param windowId id to check for
 * @return true if key is in map else false
 */
bool hasValue(std::map<int, MzSpectrum>& windowCollection, int windowId){
    return windowCollection.count(windowId);
}

/**
 * group a given spectrum into slices of fixed length along the mz axis
 * @param spec: mz spectrum to split
 * @param windowLength: length in dalton a window should have
 * @param overlapping: if true, also windows shifted by half a window length will be generated
 * @param minNumPeaks: minimum number of peaks that need to be present in a window
 * @param minIntensity: minimum intensity a peak must have to be considered
 * @return : map of grouped spectrum into windows by: windowId, mz window
 */
std::map<int, MzSpectrum> groupToWindows(const MzSpectrum& spec, double windowLength, bool overlapping=true,
                                         int minNumPeaks=1, int minIntensity=1){

    std::map<int, MzSpectrum> splits;

    for(int i = 0; i < spec.mz.size(); i++) {
        // get mz and intensity
        auto mz = spec.mz[i];
        auto intensity = spec.intensity[i];
        // calculate window key mz and intensity value should be sent to
        auto tmpKey = int(floor(mz / windowLength));

        // given key is in map
        if (hasValue(splits, tmpKey)) {
            splits[tmpKey].mz.push_back(mz);
            splits[tmpKey].intensity.push_back(intensity);
        }
            // given key is not in map
        else {
            std::vector<double> tmpMz = {mz};
            std::vector<double> tmpI = {intensity};
            splits[tmpKey] = MzSpectrum{spec.frameId, spec.scanId, tmpMz, tmpI};
        }
    }

    // calculate grouping by offset
    if(overlapping) {
        std::map<int, MzSpectrum> splitsOffset;

        for (int i = 0; i < spec.mz.size(); i++) {
            auto mz = spec.mz[i];
            auto intensity = spec.intensity[i];
            // calculate window key with offset mz and intensity value should be sent to
            auto tmpKey = -int(floor(((mz + windowLength / 2.0) / windowLength)));

            // given key is in map
            if (hasValue(splitsOffset, tmpKey)) {
                splitsOffset[tmpKey].mz.push_back(mz);
                splitsOffset[tmpKey].intensity.push_back(intensity);
            }
                // given key is not in map
            else {
                std::vector<double> tmpMz = {mz};
                std::vector<double> tmpI = {intensity};
                splitsOffset[tmpKey] = MzSpectrum{spec.frameId, spec.scanId, tmpMz, tmpI};
            }
        }
        splits.merge(splitsOffset);
    }

    std::map<int, MzSpectrum> retSplits;

    for(const auto& [bin, spectrum]: splits){
        if(spec.mz.size() >= minNumPeaks){
            auto it = max_element(std::begin(spectrum.intensity), std::end(spectrum.intensity));
            if(it >= minIntensity)
                retSplits[bin] = spectrum;
        }
    }

    return retSplits;
}

/**
 * function to bin a spectrum to a given number of decimal places
 * @param spec: spectrum to bin
 * @param resolution: number of decimal places to keep
 * @return binned: spectrum
 */
MzSpectrum binToResolution(const MzSpectrum& spec, int resolution) {

    std::map<int, int> intensityMap;
    double factor = pow(10.0, resolution);

    std::vector<double> resMz;
    std::vector<double> resI;

    for(int i = 0; i < spec.mz.size(); i++){

        // calculate binned mz value as key
        int rounded = int(roundf(spec.mz[i] * factor));
        // add intensity to generated key
        intensityMap[rounded] += spec.intensity[i];
    }

    resMz.reserve(intensityMap.size());
    resI.reserve(intensityMap.size());

    // get all mz values and sort them
    for (const auto& [key, _] : intensityMap) {
        resMz.push_back(double(key) / factor);
    }

    std::sort(resMz.begin(), resMz.end());

    // get intensity values in sorted order
    for (const auto &elem : resMz) {
        resI.push_back(intensityMap[int(elem * factor)]);
    }

    return MzSpectrum{spec.frameId, spec.scanId, resMz, resI};
}

/**
 * translate a binned mz spectrum to a sparse vector where indices are mz indices and values are intensities
 * @param spec: binned spectrum to vectorize
 * @param resolution: number of decimals given spectrum has
 * @return a sparse vectorized spectrum
 */
MzVector toVectorizedBinnedSpectrum(const MzSpectrum& spec, int resolution=2, bool normalize=false, bool sqr=false) {

    int factor = pow(10, resolution);
    std::vector<int> indices;
    indices.reserve(spec.mz.size());

    std::vector<double> values;
    values.reserve(spec.intensity.size());

    for(const auto m : spec.mz)
        indices.push_back(m * factor);

    for(const auto i : spec.intensity){
        if(normalize)
            if(sqr)
                values.push_back(sqrt(i));
            else
                values.push_back(log(i));
        else
            values.push_back(i);
    }

    return MzVector{resolution, spec.frameId, spec.scanId, indices, values};
}

#endif //TIMSLSH_SPECTRUM_H