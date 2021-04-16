// g++ -o s s.cpp -std=c++17 -Wall -W -fopenmp -O3

#include <cstdio>
#include <iostream>
#include <random>
#include <string>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <chrono>

constexpr size_t KMER_LENGTH = 5;
constexpr size_t MAX_SEQUENCES_PER_FILE = 1000000;

#include "CreateKmerSignatures.h"
#include "KMeansClustering.h"

using namespace std;
using namespace std::chrono;

const char *readLine(FILE *fp)
{
  static char lineBuffer[1048576];
  #pragma omp threadprivate(lineBuffer)
  if (fscanf(fp, "%1048575[^\n]\n", lineBuffer) < 1) {
    return nullptr;
  }
  return lineBuffer;
}

vector<string> readFastq(const char *path)
{
  vector<string> fastq;
  FILE *fp = fopen(path, "r");
  if (!fp) {
    fprintf(stderr, "Failed to open %s\n", path);
    exit(1);
  }
  for (size_t i = 0; i < MAX_SEQUENCES_PER_FILE; i++) {
    if (!readLine(fp)) break;
    const char *seq = readLine(fp);
    fastq.push_back(string(seq));
    if (!readLine(fp)) break;
    if (!readLine(fp)) break;
  }
  fclose(fp);
  return fastq;
}

vector<vector<double>> createDocumentSignatures(const vector<string> &fastq, const KmerSignatures &kmerSignatures)
{
  fprintf(stderr, "Creating signatures for %zu sequences (expected size: %zumb)\n", fastq.size(), fastq.size() * kmerSignatures.getDimensionality() * sizeof(double) / 1024 / 1024);
  vector<vector<double>> docSigs(fastq.size(), vector<double>(kmerSignatures.getDimensionality(), 0.0));
  #pragma omp parallel for
  for (size_t i = 0; i < fastq.size(); i++) {
    kmerSignatures.getDocumentSignature(fastq[i], docSigs[i]);
  }
  
  return docSigs;
}

//                                1      2      3      4      5      6      7      8
const vector<string> CITIES = {"ARN", "BCN", "BER", "DEN", "DOH", "FAI", "HKG", "ICN",
// 9     10     11     12     13     14     15     16     17     18     19     20     21     22     23
"IEV", "ILR", "KUL", "LCY", "LIS", "NYC", "OFF", "SAO", "SCL", "SDJ", "SFO", "SGP", "TPE", "TYO", "ZRH"};

int getCityNum(const char *path)
{
  for (size_t i = 0; i < CITIES.size(); i++) {
    string pattern = "_" + CITIES[i] + "_";
    if (strstr(path, pattern.c_str())) return i + 1;
  }
  return 0;
}

string filename(const char *path)
{
  const char *lastSlash = strrchr(path, '/');
  if (!lastSlash) return string(path);
  return string(lastSlash + 1);
}

int main(int argc, char **argv)
{
  vector<vector<string>> fastqs(argc);
  
  int c = 0;
  #pragma omp parallel for num_threads(8)
  for (int arg = 1; arg < argc; arg++) {
    fastqs[arg] = readFastq(argv[arg]);
    #pragma omp critical
    {
      c++;
      fprintf(stderr, "%s %d/%d\n", argv[arg], c, argc-1); 
    }
  }
  
  bool first = true;
  for (;;) {
    duration<double> total_sig_time(0);
    duration<double> total_clustering_time(0);
    auto time_before = steady_clock::now();
    char prefix[100];
    size_t kmerLength, clusterCount, desiredSigSize;
    if (scanf("%s %zu %zu %zu\n", prefix, &kmerLength, &clusterCount, &desiredSigSize) < 4) {
      if (first) {
        fprintf(stderr, "Need to pass prefix, kmer length and cluster count via stdin\n");
        exit(1);
      } else {
        break;
      }
    }
    first = false;
    
    KmerSignatures kmerSignatures(kmerLength, desiredSigSize);
    
    vector<int> sampleNums;
    vector<int> cityNums;
    vector<double> clusterMeans;
    vector<string> samplePaths;
    samplePaths.push_back("(no sample)");
    
    for (int arg = 1; arg < argc; arg++) {
      samplePaths.push_back(filename(argv[arg]));
    }
    
    //#pragma omp parallel for schedule(dynamic)
    for (int arg = 1; arg < argc; arg++) {
      const char *path = argv[arg];
      const auto &fastq = fastqs[arg];
      int cityNum = getCityNum(path);
      fprintf(stderr, "%s - city num: %d\n", path, cityNum);
      fprintf(stderr, "Read %zu sequences\n", fastq.size());
      auto time_before_sigcreation = steady_clock::now();
      auto sigs = createDocumentSignatures(fastq, kmerSignatures);
      auto time_after_sigcreation = steady_clock::now();
      auto time_taken_sigcreation = duration_cast<duration<double>>(time_after_sigcreation - time_before_sigcreation);
      total_sig_time += time_taken_sigcreation;
      
      fprintf(stderr, "Created %zu signatures\n", sigs.size());
      
      vector<size_t> clusters;
      vector<vector<double>> outputClusterMeans;
      auto time_before_clustering = steady_clock::now();
      clusterSignaturesX(sigs, clusters, outputClusterMeans, clusterCount, 100);
      auto time_after_clustering = steady_clock::now();
      auto time_taken_clustering = duration_cast<duration<double>>(time_after_clustering - time_before_clustering);
      total_clustering_time += time_taken_clustering;
      fprintf(stderr, "Created %zu clusters\n", outputClusterMeans.size());
      
      //#pragma omp critical
      for (size_t i = 0; i < outputClusterMeans.size(); i++) {
        sampleNums.push_back(arg);
        cityNums.push_back(cityNum);
        for (double d : outputClusterMeans[i]) {
          clusterMeans.push_back(d);
        }
      }
    }
    
    auto time_after = steady_clock::now();
    auto time_taken = duration_cast<duration<double>>(time_after - time_before);
    
    size_t signatureWidth = kmerSignatures.getDimensionality();
    size_t numClusters = sampleNums.size();
    
    char outputPath[1024];
    sprintf(outputPath, "out/%s%zu_k%zu_%zu.time", prefix, clusterCount, kmerLength, signatureWidth);
    FILE *fo = fopen(outputPath, "w");
    fprintf(fo, "Total time: %f seconds\n", time_taken.count());
    fprintf(fo, "Signature creation time: %f seconds\n", total_sig_time.count());
    fprintf(fo, "Clustering time: %f seconds\n", total_clustering_time.count());
    fclose(fo);
    sprintf(outputPath, "out/%s%zu_k%zu_%zu.sigs", prefix, clusterCount, kmerLength, signatureWidth);
    fo = fopen(outputPath, "wb");
    
    fwrite(&signatureWidth, sizeof(size_t), 1, fo);
    fwrite(&numClusters, sizeof(size_t), 1, fo);
    fwrite(sampleNums.data(), sizeof(int), sampleNums.size(), fo);
    fwrite(cityNums.data(), sizeof(int), cityNums.size(), fo);
    fwrite(clusterMeans.data(), sizeof(double), clusterMeans.size(), fo);
    
    for (string path : samplePaths) {
      fprintf(fo, "%s\n", path.c_str());
    }
    fclose(fo);
  }
}
