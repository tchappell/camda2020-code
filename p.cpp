// g++ -o p p.cpp -std=c++11 -fopenmp -Wall -O3

#include <cstdio>
#include <cstring>
#include <iostream>
#include <random>
#include <string>
#include <map>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <vector>
#include <chrono>

using namespace std;
using namespace std::chrono;

double getDist(const double *a, const double *b, size_t dim)
{
  double dist = 0.0;
  for (size_t i = 0; i < dim; i++) {
    double d = a[i] - b[i];
    dist += d * d;
  }
  
  return dist;
}

void fmterr() {
  fprintf(stderr, "Format error\n");
  exit(1);
}

//                                     1      2      3      4      5      6      7      8
const vector<string> CITIES = {"???", "ARN", "BCN", "BER", "DEN", "DOH", "FAI", "HKG", "ICN",
// 9     10     11     12     13     14     15     16     17     18     19     20     21     22     23
"IEV", "ILR", "KUL", "LCY", "LIS", "NYC", "OFF", "SAO", "SCL", "SDJ", "SFO", "SGP", "TPE", "TYO", "ZRH"};
const vector<pair<string, int>> MYSTERY = {{"VxRfsZ5CEipk", 22}, {"skcahSj00d0n", 22}, {"yVuBUyFA3BaG", 7}, {"CAov7ffceNbk", 7}, {"rFns8qLWHkEf", 7}, {"hB7CJ4Jbvo9o", 7}, {"5wCWl46pJ3MX", 21}, {"VJcXyLAgaaSw", 21}, {"EA3qOpmxgoy9", 22}, {"skgN1DPY3tfx", 7}, {"QMsho8YSU54N", 7}, {"mUWBc0gzBK0r", 22}, {"PEIQYIGgdqRR", 22}, {"7mm7F8ebohvI", 7}, {"a0pTYjAUGvLb", 7}, {"bhDkIhhbsn63", 7}, {"7XChw0bYwSvi", 7}, {"hiBLpdVbSkOh", 7}, {"0V1JrUr3qPfY", 7}, {"vRfWohxWD0Hk", 7}, {"klKodXilNzwr", 7}, {"hYlkiU9h5dAk", 7}, {"ONPzTKdGTXuB", 21}, {"fSsPSaDvTDFo", 21}, {"my9YP1uGBLAy", 21}, {"16nKqloPqYql", 21}, {"B9TwLHfoYhG2", 21}, {"w6s878eMcnQa", 21}, {"mxb4uDYFG1bv", 21}, {"jYCnTTvl48gm", 21}, {"xLuyUbkUQnQD", 21}, {"3CtKjudk0Wt2", 22}, {"TCtAKW6dGptn", 22}, {"d0z9i0tIWKPV", 22}, {"0Ae1dHGe0DyY", 22}, {"TmoWLQAjWb1V", 22}, {"MNjYO4DHYxuY", 22}, {"PaJxnodc0WTC", 22}, {"D92WwTPqqCN9", 22}, {"uEFSltwT1KjW", 22}, {"PRIQj2M2oiJQ", 23}, {"gDu37tpc22jp", 23}, {"ThvUPMBYvXVR", 23}, {"aJvQ1zlMspQg", 23}, {"Tbw0MUSK9HWs", 23}, {"o7RBip4p0ah1", 23}, {"uhP8ZXxgiWms", 23}, {"rwVA6i8hcmVZ", 23}, {"ANCJV9G3edsg", 23}, {"hS8PZiES2QLb", 23}, {"j9lzr0Wx8GUE", 23}, {"5lEou6P7Wfid", 23}, {"Enr7963z6VQq", 23}, {"a0ghpkzgxC4G", 23}, {"6YwsuU6EHlT0", 9}, {"Begi5ms9HQE6", 9}, {"gfjfi52jRl7W", 9}, {"N6HxsNCeOpCu", 9}, {"O3Byrj7fOgpd", 9}, {"r2nmgQkfTTVo", 9}, {"RXA8KOX22npi", 9}, {"sUMZQLn3yBSB", 9}, {"syatOMVYSSCU", 9}, {"VoFWzgsqkR2W", 9}, {"Vvw1JWzmLMWr", 9}
};

int mysteryLookup(const char *path)
{
  for (const auto &m : MYSTERY) {
    if (strstr(path, m.first.c_str())) {
      return m.second;
    }
  }
  return 0;
}

struct SignatureFile {
  size_t signatureWidth;
  size_t numClusters;
  vector<int> sampleNums;
  vector<int> cityNums;
  vector<double> clusterMeans;
  vector<string> samplePaths;
  SignatureFile(const char *path) {
    FILE *fp = fopen(path, "rb");
    if (!fp) {
      fprintf(stderr, "Failed to open %s\n", path);
      exit(1);
    }
    
    if (fread(&signatureWidth, sizeof(size_t), 1, fp) != 1) fmterr();
    if (fread(&numClusters, sizeof(size_t), 1, fp) != 1) fmterr();  
    sampleNums.resize(numClusters);
    cityNums.resize(numClusters);
    clusterMeans.resize(numClusters * signatureWidth);
    if (fread(sampleNums.data(), sizeof(int), sampleNums.size(), fp) != numClusters) fmterr();
    if (fread(cityNums.data(), sizeof(int), cityNums.size(), fp) != numClusters) fmterr();
    if (fread(clusterMeans.data(), sizeof(double), clusterMeans.size(), fp) != numClusters * signatureWidth) fmterr();

    samplePaths.clear();
    
    for (;;) {
      char buf[1024];
      if (fscanf(fp, "%[^\n]\n", buf) < 1) break;
      samplePaths.push_back(string(buf));
    }
    
    for (size_t i = 0; i < samplePaths.size(); i++) {
      auto sample = samplePaths[i];
      int m = mysteryLookup(sample.c_str());
      if (m != 0) {
        //fprintf(stderr, "%s -> %s\n", sample.c_str(), CITIES[m].c_str());
        // Now we know that every cluster with a sampleNum of i is from city m
        for (size_t j = 0; j < numClusters; j++) {
          if (sampleNums[j] == static_cast<int>(i)) {
            cityNums[j] = m;
          }
        }
      }
    }
    fclose(fp);
  }
};


void doEval(const char *path, SignatureFile &train, SignatureFile &test, bool verbose)
{
    auto time_before_search = steady_clock::now();

    vector<int> predictedSample(test.numClusters);
    bool leaveOneOut = false;
    if (&train == &test) {
      leaveOneOut = true;
    }
    if (train.signatureWidth != test.signatureWidth) {
      fprintf(stderr, "Mismatched signature sizes (%zu != %zu)\n", train.signatureWidth, test.signatureWidth);
      exit(1);
    }
    
    #pragma omp parallel for
    for (size_t i = 0; i < test.numClusters; i++) {
      // Do a nearest neighbour search for this (test) cluster, but only within this sample
      int sampleNum = test.sampleNums[i];
      size_t closestCluster = numeric_limits<size_t>::max();
      double closestDist = numeric_limits<double>::max();
      double *source_signature = &test.clusterMeans[i * test.signatureWidth];
      
      for (size_t j = 0; j < train.numClusters; j++) {
        if (leaveOneOut && train.sampleNums[j] == sampleNum) continue;
        double *dest_signature = &train.clusterMeans[j * train.signatureWidth];
        double d = getDist(source_signature, dest_signature, train.signatureWidth);
        if (d < closestDist) {
          closestDist = d;
          closestCluster = j;
        }
      }
      predictedSample[i] = train.cityNums[closestCluster];
    }
    
    map<int, unordered_map<int, int>> predictions;
    unordered_map<int, int> actualCity;
    
    for (size_t i = 0; i < test.numClusters; i++) {
      int sampleNum = test.sampleNums[i];
      int predictedCity = predictedSample[i];
      predictions[sampleNum][predictedCity]++;
      actualCity[sampleNum] = test.cityNums[i];
    }
    
    
    int predictionCounts = 0;
    int correct = 0;
    int samplePathCount = test.samplePaths.size();
    for (const auto &sample_prediction : predictions) {
      int sampleNum = sample_prediction.first;
      int correctCity = actualCity[sampleNum];
      
      int cityMax = 0;
      int predictedCity = 0;
      int predictedCityCount = 0;
      for (const auto &prediction : sample_prediction.second) {
        int predCounts = prediction.second;
        int predCity = prediction.first;
        //fprintf(stderr, "%s %d\n", CITIES[predCity].c_str(), predCounts);

        if (predCounts > cityMax) {
          cityMax = predCounts;
          predictedCity = predCity;
          predictedCityCount = 0;
        }
        if (predCounts == cityMax) {
          predictedCityCount++;
        }
      }
      if (correctCity != 0) {
        predictionCounts++;
        if (sample_prediction.second.count(correctCity) == 1 && sample_prediction.second.at(correctCity) == cityMax) {
          if (predictedCityCount == 1) {
            //if (predictedCity == correctCity) {
            correct++;
          } else {
            //fprintf(stderr, "Correct but thrown out due to predictedCityCount = %d\n", predictedCityCount);
          }
        }
      }
      
      if (verbose) {
        if (sampleNum < samplePathCount) {
          fprintf(stderr, "%s", test.samplePaths[sampleNum].c_str());
        } else {
          fprintf(stderr, "Sample %d", sampleNum);
        }
        fprintf(stderr, " - %s (%s)\n", CITIES[predictedCity].c_str(), CITIES[correctCity].c_str());
      }
    }
    auto time_after_search = steady_clock::now();
    auto time_taken_search = duration_cast<duration<double>>(time_after_search - time_before_search);

    double acc = 100.0 * correct / predictionCounts;
    fprintf(stderr, "%s: Correctly predicted %d/%d (%f%%) (%f seconds)\n", path, correct, predictionCounts, acc, time_taken_search.count());

}

int evaluate(int argc, char **argv) 
{
  for (int arg = 2; arg < argc; arg += 2) {
    auto trainSignature = SignatureFile(argv[arg]);
    auto testSignature = SignatureFile(argv[arg + 1]);
    
    doEval(argv[arg], trainSignature, testSignature, false);
  }
  return 0;
}

int main(int argc, char **argv)
{
  if (argc < 2) {
    fprintf(stderr, "Usage: %s (-e) [sigs file...]\n", argv[0]);
    fprintf(stderr, "If -e is used, this program will expect pairs of sigs files to follow\n");
    exit(1);
  }
  if (strcmp(argv[1], "-e")==0) {
    return evaluate(argc, argv);
  }
  
  for (int arg = 1; arg < argc; arg++) {
    auto signature = SignatureFile(argv[arg]);
    
    doEval(argv[arg], signature, signature, false);
  }
  return 0;
}
