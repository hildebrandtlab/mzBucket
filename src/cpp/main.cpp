#include <set>
#include <iostream>
#include <map>
#include <vector>
#include <cmath>
#include <functional>
#include <string>
#include <fstream>
#include "cxxopts.hpp"
#include "spectrum.h"
#include "hashing.h"
#include </home/dateschn/opentims/opentims++/opentims_all.cpp>
//#include </Users/davidteschner/Documents/promotion/opentims/opentims++/opentims_all.cpp>

int main(int argc, char *argv[]) {

    cxxopts::Options options("mzBucket", "Fast Cosim-Hashes for TimsTOF Spectra");
    options.add_options()
        ("f,frameId", "frame number", cxxopts::value<int>()->default_value("1"))
            ("w, windowlength", "length of subspectra(da)", cxxopts::value<double>()->default_value("10"))
            ("o, overlapping", "whether windows overlap", cxxopts::value<std::string>()->default_value("true"))
            ("k, OR", "number of ORs", cxxopts::value<int>())
            ("l, AND", "number of ANDs", cxxopts::value<int>())
            ("c, csvPath", "a csv file for testing", cxxopts::value<std::string>())
            ("d, brukerDataSet", "path to bruker dataset", cxxopts::value<std::string>())
            ("b, binaryPath", "path to bruker binary", cxxopts::value<std::string>())
            ("n, normalize", "if true, norm(I) will be used", cxxopts::value<std::string>()->default_value("false"))
            ("s, sqrt", "if true, sqrt norm will be used", cxxopts::value<std::string>()->default_value("false"))
            ("t, threads", "number of threads", cxxopts::value<int>()->default_value("4"))
            ("h, help", "Print usage");

    auto result = options.parse(argc, argv);

    if (result.count("help"))
    {
        std::cout << options.help() << std::endl;
        exit(0);
    }

    if(!result.count("frameId") || !result.count("AND") || !result.count("OR")){
        std::cout << options.help() << std::endl;
        exit(0);
    }

    if(result.count("csvPath")  && result.count("brukerDataSetPath") && result.count("brukerBinaryPath")){
        std::cout << options.help() << std::endl;
        exit(0);
    }

    int frameId = result["frameId"].as<int>();
    double windowlength = result["windowlength"].as<double>();
    std::string overlappingS = result["overlapping"].as<std::string>();
    std::string normalizeS = result["normalize"].as<std::string>();
    std::string sqrtS = result["sqrt"].as<std::string>();
    int k = result["AND"].as<int>();
    int l = result["OR"].as<int>();
    int numThreads = result["threads"].as<int>();

    bool overlapping = false;
    bool normalize = false;
    bool sqrt = false;

    if(overlappingS == "true")
        overlapping = true;

    if(normalizeS == "true")
        normalize = true;

    if(sqrtS == "true")
        sqrt = true;

    bool isTestCase = false;

    if(result.count("csvPath"))
        isTestCase = true;

    std::cout << "LSH SETTINGS: " << std::endl;
    std::cout << "______________________________" << std::endl;
    std::cout << "window length(dalton): " << windowlength << std::endl;
    std::cout << "windows overlap      : " << overlappingS << std::endl;
    std::cout << "number of ANDs       : " << k << std::endl;
    std::cout << "number of ORs        : " << l << std::endl;
    std::cout << "normalize intensities: " << normalizeS << std::endl;

    if (normalize)
        if(sqrt)
            std::cout << "normalize function   : " << "SQRT" << std::endl;
        else
            std::cout << "normalize function   : " << "LOG" << std::endl;

    std::cout << "______________________________" << std::endl;
    std::cout << "" << std::endl;

    // test mode read data from csv
    if(isTestCase){

        auto csvPath = result["csvPath"].as<std::string>();

        // file object for reading
        std::ifstream file(csvPath);

        // skip the first line
        std::string line;
        std::getline(file, line);

        std::stringstream stringstream;
        stringstream.str(line);
        std::string substring;
        std::vector<int> scans;
        std::vector<double> mzs;
        std::vector<double> intensitys;
        std::vector<std::string> labels;
        std::vector<std::string> row;

        while (std::getline(file, line)) {
            substring.clear();
            stringstream.clear();
            stringstream.str(line);
            while (std::getline(stringstream, substring, ',')) {
                row.push_back(substring);
            }

            auto sca = std::atoi(row[0].c_str());
            auto m = std::atof(row[1].c_str());
            auto intens = std::atof(row[2].c_str());
            auto labe = row[3];

            scans.push_back(sca);
            mzs.push_back(m);
            intensitys.push_back(intens);
            labels.push_back(labe);

            row.clear();
        }

        file.close();

        std::map<int, MzSpectrum> specMap;
        std::vector<std::string> la;

        for(auto peak_id = 0; peak_id < scans.size(); peak_id++){
            auto scan = scans[peak_id];
            auto mz = mzs[peak_id];
            auto i = intensitys[peak_id];
            if(hasValue(specMap, scan)){
                specMap[scan].mz.push_back(mz);
                specMap[scan].intensity.push_back(i);
            }
            else {
                specMap[scan] = MzSpectrum{1, static_cast<int>(scan), {mz}, {i}};
            }
        }

        auto checkSet = getHashes(specMap, numThreads, windowlength, overlapping, k, l, normalize, sqrt);

        for(auto peak_id = 0; peak_id < mzs.size(); peak_id++){

           auto scan = scans[peak_id];
           auto mz = mzs[peak_id];
           auto i = intensitys[peak_id];

           auto bin = int(floor(mzs[peak_id] / windowlength));
           auto oeBin = -int(floor(((mzs[peak_id] + windowlength / 2.0) / windowlength)));

           if(hasValue(checkSet, scan)){
               auto s = checkSet[scan];
               if(s.count(bin) > 0 || s.count(oeBin) > 0){
                   la.push_back(labels[peak_id]);
               }
           }
        }

       auto tCounter = 0;
       auto fCounter = 0;

       auto tCounterReal = 0;
       auto fCounterReal = 0;

       for(const auto s: labels){
           if(s == "True")
               tCounterReal += 1;
           else
               fCounterReal += 1;
       }

       for(const auto s: la){
           if(s == "True")
               tCounter += 1;
           else
               fCounter += 1;
       }

       std::cout << "Real True: " << tCounterReal << " | Real False: " << fCounterReal << std::endl;
       std::cout << "Pred True: " << tCounter << " | Pred False: " << fCounter << std::endl;
   }
   else{

       auto brukerDataSetPath = result["brukerDataSet"].as<std::string>();
       auto brukerBinaryPath = result["binaryPath"].as<std::string>();

       DefaultTof2MzConverterFactory::setAsDefault<BrukerTof2MzConverterFactory, const char*>(brukerBinaryPath.c_str());
       DefaultScan2InvIonMobilityConverterFactory::setAsDefault<BrukerScan2InvIonMobilityConverterFactory, const char*>(brukerBinaryPath.c_str());

       // Open the dataset
       TimsDataHandle TDH(brukerDataSetPath);

       // Allocate buffers for data: instead of reallocating for every frame, we just allocate a buffer that will fit all
       // and reuse it.
       const size_t buffer_size_needed = TDH.max_peaks_in_frame();

       std::unique_ptr<uint32_t[]> scan_ids = std::make_unique<uint32_t[]>(buffer_size_needed);
       std::unique_ptr<uint32_t[]> intensities = std::make_unique<uint32_t[]>(buffer_size_needed);

       std::unique_ptr<double[]> mzs = std::make_unique<double[]>(buffer_size_needed);

       std::unique_ptr<double[]> retention_times = std::make_unique<double[]>(buffer_size_needed);

       TimsFrame& frame = TDH.get_frame(frameId);
       frame.save_to_buffs(nullptr, scan_ids.get(), nullptr,
                           intensities.get(), mzs.get(), nullptr, nullptr);


       std::map<int, MzSpectrum> specMap;

       for(size_t peak_id = 0; peak_id < frame.num_peaks; peak_id++){
           auto scan = scan_ids[peak_id];
           auto mz = mzs[peak_id];
           auto i = intensities[peak_id];
           if(hasValue(specMap, scan)){
               specMap[scan].mz.push_back(mz);
               specMap[scan].intensity.push_back(i);
           }
           else {
               specMap[scan] = MzSpectrum{frameId, static_cast<int>(scan), {mz}, {double(i)}};
           }
       }

       auto checkSet = getHashes(specMap, numThreads, windowlength, overlapping, k, l, normalize, sqrt);

       std::vector<double> mzValues;
       std::vector<int> intens;
       std::vector<int> scans;

       std::vector<std::string> la;

           for(size_t peak_id = 0; peak_id < frame.num_peaks; peak_id++){

               auto scan = scan_ids[peak_id];
               auto mz = mzs[peak_id];
               auto i = mzs[peak_id];

               auto bin = int(floor(mzs[peak_id] / windowlength));
               auto oeBin = -int(floor(((mzs[peak_id] + windowlength / 2.0) / windowlength)));

               if(hasValue(checkSet, scan)){
                   auto s = checkSet[scan];
                   if(s.count(bin) > 0 || s.count(oeBin) > 0){
                       mzValues.push_back(mz);
                       intens.push_back(i);
                       scans.push_back(scan);
                   }
               }
           }

           std::cout << "Number peaks raw     : " << frame.num_peaks << std::endl;
           std::cout << "Number Peaks filtered: " << mzValues.size() << std::endl;
   }
}
