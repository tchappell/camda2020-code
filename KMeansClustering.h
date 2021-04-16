#include <cstdio>
#include <cstdint>
#include <iostream>
#include <random>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <algorithm>

using namespace std;

template<class RNG>
vector<vector<double>> createRandomSigs(RNG &&rng, const vector<vector<double>> &sigs, size_t clusterCount)
{
  vector<vector<double>> clusterSigs;
  
  uniform_int_distribution<size_t> dist(0, sigs.size() - 1);
  set<vector<double>> uniqueSigs;
  
  while (uniqueSigs.size() < clusterCount) {
    size_t sigIdx = dist(rng);
    uniqueSigs.insert(sigs[sigIdx]);
  }

  for (const auto &sig : uniqueSigs) {
    clusterSigs.push_back(sig);
  }
  
  return clusterSigs;
}

template<class RNG>
vector<uint64_t> createRandomSigs(size_t signatureSize, RNG &&rng, const vector<uint64_t> &sigs, size_t clusterCount)
{
  vector<uint64_t> clusterSigs(signatureSize * clusterCount);
  size_t signatureCount = sigs.size() / signatureSize;
  uniform_int_distribution<size_t> dist(0, signatureCount - 1);
  bool finished = false;
  
  unordered_set<string> uniqueSigs;
  for (size_t i = 0; i < signatureCount; i++) {
    size_t sig = dist(rng);
    string sigData(signatureSize * sizeof(uint64_t), ' ');
    memcpy(&sigData[0], &sigs[sig * signatureSize], signatureSize * sizeof(uint64_t));
    uniqueSigs.insert(sigData);
    if (uniqueSigs.size() >= clusterCount) {
      finished = true;
      break;
    }
  }
  
  if (!finished) {
    fprintf(stderr, "Error: there appear to be fewer unique signatures than clusters (%zu uniques).\n", uniqueSigs.size());
    fprintf(stderr, "Reduce the cluster count\n");\
    
    vector<uint64_t> binSig(signatureSize);
    for (const auto &sig : uniqueSigs) {
      memcpy(binSig.data(), sig.data(), signatureSize * sizeof(uint64_t));
      for (size_t i = 0; i < signatureSize; i++) {
        printf("%0lx", binSig[i]);
      }
      printf("\n");
    }
    exit(1);
  }
  
  size_t i = 0;
  for (const auto &sig : uniqueSigs) {
    memcpy(&clusterSigs[i * signatureSize], sig.data(), signatureSize * sizeof(uint64_t));
    i++;
  }
  
  if (i != uniqueSigs.size()) {
    fprintf(stderr, "Mismatch: %zu != %zu\n", i, uniqueSigs.size());
    exit(1);
  }
  
  return clusterSigs;
}

double reclusterSignatures(vector<size_t> &clusters, const vector<vector<double>> &meanSigs, const vector<vector<double>> &sigs, size_t clusterCount)
{
  vector<double> hdist(clusters.size());
  #pragma omp parallel for
  for (size_t sig = 0; sig < clusters.size(); sig++) {
    const auto sourceSignature = sigs[sig];
    size_t minHdCluster = 0;
    double minHd = numeric_limits<double>::max();
    for (size_t cluster = 0; cluster < clusterCount; cluster++) {
      const auto &clusterSignature = meanSigs[cluster];
      double hd = 0.0;
      for (size_t i = 0; i < sourceSignature.size(); i++) {
        double d = sourceSignature[i] - clusterSignature[i];
        hd += d * d;
      }
      if (hd < minHd) {
        minHd = hd;
        minHdCluster = cluster;
      }
    }
    hdist[sig] = minHd;
    clusters[sig] = minHdCluster;
  }
  
  double mse = 0.0;
  for (size_t sig = 0; sig < clusters.size(); sig++) {
    mse += hdist[sig] * hdist[sig];
  }
  return mse;
}

size_t reclusterSignatures(size_t signatureSize, vector<size_t> &clusters, const vector<uint64_t> &meanSigs, const vector<uint64_t> &sigs, size_t clusterCount)
{
  vector<size_t> hdist(clusters.size());
  #pragma omp parallel for
  for (size_t sig = 0; sig < clusters.size(); sig++) {
    const uint64_t *sourceSignature = &sigs[sig * signatureSize];
    size_t minHdCluster = 0;
    size_t minHd = numeric_limits<size_t>::max();
    for (size_t cluster = 0; cluster < clusterCount; cluster++) {
      const uint64_t *clusterSignature = &meanSigs[cluster * signatureSize];
      size_t hd = 0;
      for (size_t i = 0; i < signatureSize; i++) {
        hd += __builtin_popcountll(sourceSignature[i] ^ clusterSignature[i]);
      }
      if (hd < minHd) {
        minHd = hd;
        minHdCluster = cluster;
      }
    }
    hdist[sig] = minHd;
    clusters[sig] = minHdCluster;
  }
  
  size_t mse = 0;
  for (size_t sig = 0; sig < clusters.size(); sig++) {
    mse += hdist[sig] * hdist[sig];
  }
  return mse;
}

vector<vector<double>> createClusterSigs(const vector<vector<size_t>> &clusterLists, const vector<vector<double>> &sigs)
{
  size_t signatureWidth = sigs[0].size();
  vector<vector<double>> clusterSigs(clusterLists.size());
  #pragma omp parallel
  {
    vector<double> unflattenedSignature(signatureWidth);
    #pragma omp for
    for (size_t cluster = 0; cluster < clusterLists.size(); cluster++) {
      fill(begin(unflattenedSignature), end(unflattenedSignature), 0.0);
      
      for (size_t signature : clusterLists[cluster]) {
        const auto signatureData = sigs[signature];
        for (size_t i = 0; i < signatureWidth; i++) {
          unflattenedSignature[i] += signatureData[i];
        }
      }
      
      for (size_t i = 0; i < signatureWidth; i++) {
        unflattenedSignature[i] /= clusterLists[cluster].size();
      }
      clusterSigs[cluster] = unflattenedSignature;
    }
  }
  return clusterSigs;
}

vector<uint64_t> createClusterSigs(size_t signatureSize, const vector<vector<size_t>> &clusterLists, const vector<uint64_t> &sigs)
{
  size_t clusterCount = clusterLists.size();
  size_t signatureWidth = signatureSize * 64;
  vector<uint64_t> clusterSigs(signatureSize * clusterCount);
  #pragma omp parallel
  {
    vector<int> unflattenedSignature(signatureWidth);
    #pragma omp for
    for (size_t cluster = 0; cluster < clusterLists.size(); cluster++) {
      fill(begin(unflattenedSignature), end(unflattenedSignature), 0);
      
      for (size_t signature : clusterLists[cluster]) {
        const uint64_t *signatureData = &sigs[signatureSize * signature];
        for (size_t i = 0; i < signatureWidth; i++) {
          uint64_t signatureMask = (uint64_t)1 << (i % 64);
          if (signatureMask & signatureData[i / 64]) {
            unflattenedSignature[i] += 1;
          } else {
            unflattenedSignature[i] -= 1;
          }
        }
      }
      
      uint64_t *flattenedSignature = &clusterSigs[cluster * signatureSize];
      for (size_t i = 0; i < signatureWidth; i++) {
        if (unflattenedSignature[i] > 0) {
          flattenedSignature[i / 64] |= (uint64_t)1 << (i % 64);
        }
      }
    }
  }
  return clusterSigs;
}

vector<vector<size_t>> createClusterLists(const vector<size_t> &clusters, size_t clusterCount)
{
  vector<vector<size_t>> clusterLists(clusterCount);
  for (size_t i = 0; i < clusters.size(); i++) {
    clusterLists[clusters[i]].push_back(i);
  }
  return clusterLists;
}

void clusterSignatures(const vector<vector<double>> &sigs, vector<size_t> &clusters, vector<vector<double>> &outputClusterMeans, size_t clusterCount, int kMeansIterations)
{
  auto rng = ranlux24_base();
  
  clusters.resize(sigs.size());
  vector<vector<size_t>> clusterLists;
  vector<vector<double>> meanSigs = createRandomSigs(rng, sigs, clusterCount);
  double lastMse = numeric_limits<double>::max();
  for (int iteration = 0; iteration < kMeansIterations; iteration++) {
    fprintf(stderr, "Iteration %d", iteration);
    double mse = reclusterSignatures(clusters, meanSigs, sigs, clusterCount);
    fprintf(stderr, "  (mse: %f)\n", mse);
    clusterLists = createClusterLists(clusters, clusterCount);
    
    meanSigs = createClusterSigs(clusterLists, sigs);
    if (mse * 1.000001 > lastMse) {
      fprintf(stderr, "Stopped converging\n");
      break;
    }
    lastMse = mse;
  }
  
  size_t nonEmptyClusters = 0;
  for (const vector<size_t> &clusterList : clusterLists) {
    if (!clusterList.empty()) {
      nonEmptyClusters++;
    }
  }
  
  outputClusterMeans = move(meanSigs);
  fprintf(stderr, "%llu/%llu non-empty clusters\n", static_cast<unsigned long long>(nonEmptyClusters), static_cast<unsigned long long>(clusterCount));
}

void clusterSignatures(const vector<uint64_t> &sigs, size_t sigSize, vector<size_t> &clusters, vector<uint64_t> &outputClusterMeans, size_t clusterCount, int kMeansIterations)
{
  size_t sigCount = sigs.size() / sigSize;
  
  auto rng = ranlux24_base();
  
  clusters.resize(sigCount);
  vector<vector<size_t>> clusterLists;
  vector<uint64_t> meanSigs = createRandomSigs(sigSize, rng, sigs, clusterCount);
  size_t lastMse = numeric_limits<size_t>::max();
  for (int iteration = 0; iteration < kMeansIterations; iteration++) {
    fprintf(stderr, "Iteration %d", iteration);
    size_t mse = reclusterSignatures(sigSize, clusters, meanSigs, sigs, clusterCount);
    fprintf(stderr, "  (mse: %zu)\n", mse);
    clusterLists = createClusterLists(clusters, clusterCount);
    
    meanSigs = createClusterSigs(sigSize, clusterLists, sigs);
    if (mse >= lastMse) {
      fprintf(stderr, "Stopped converging\n");
      break;
    }
    lastMse = mse;
  }
  
  size_t nonEmptyClusters = 0;
  for (const vector<size_t> &clusterList : clusterLists) {
    if (!clusterList.empty()) {
      nonEmptyClusters++;
    }
  }
  
  outputClusterMeans = move(meanSigs);
  fprintf(stderr, "%llu/%llu non-empty clusters\n", static_cast<unsigned long long>(nonEmptyClusters), static_cast<unsigned long long>(clusterCount));
}

void tsvqSignatures(const vector<vector<double>> &sigs, vector<size_t> &clusters, vector<vector<double>> &outputClusterMeans, size_t clusterCount, int kMeansIterations)
{
  if (clusterCount != 100) {
    fprintf(stderr, "Error: tsvqSignatures only works with c=100\n");
    exit(1);
  }
  // Hardcoded to handle 100 clusters. We will do 10 then 10.
  vector<vector<double>> clMeans;
  vector<size_t> initialClusters;

  clusterSignatures(sigs, initialClusters, clMeans, 10, kMeansIterations);
  
  clusters.resize(sigs.size());
  outputClusterMeans.resize(100);
  for (size_t i = 0; i < 10; i++) {
    vector<size_t> originalPos;
    vector<vector<double>> local_sigs;
    
    for (size_t j = 0; j < sigs.size(); j++) {
      if (initialClusters[j] == i) {
        local_sigs.push_back(move(sigs[j]));
        originalPos.push_back(j);
      }
    }
    
    vector<vector<double>> local_means;
    vector<size_t> local_clusters;
    clusterSignatures(local_sigs, local_clusters, local_means, 10, kMeansIterations);
    for (size_t j = 0; j < local_means.size(); j++) {
      outputClusterMeans[i * 10 + j] = move(local_means[j]);
    }
    for (size_t j = 0; j < local_sigs.size(); j++) {
      clusters[originalPos[j]] = local_clusters[j];
    }
  }
}

void clusterToOne(const vector<vector<double>> &sigs, vector<size_t> &clusters, vector<vector<double>> &outputClusterMeans)
{
  size_t sigSize = sigs[0].size();
  vector<double> meanSig(sigSize, 0.0);
  #pragma omp parallel
  {
    vector<double> local_totalSig(sigSize, 0.0);
    #pragma omp for
    for (size_t i = 0; i < sigs.size(); i++) {
      const auto &sig = sigs[i];
      for (size_t j = 0; j < sigSize; j++) {
        local_totalSig[j] += sig[j];
      }
    }
    #pragma omp critical
    {
      for (size_t j = 0; j < sigSize; j++) {
        meanSig[j] += local_totalSig[j];
      }
    }
  }
  for (size_t j = 0; j < sigSize; j++) {
    meanSig[j] /= sigs.size();
  }
  clusters = vector<size_t>(sigs.size(), 0);
  outputClusterMeans.clear();
  outputClusterMeans.push_back(meanSig);
}

void clusterSignaturesX(const vector<vector<double>> &sigs, vector<size_t> &clusters, vector<vector<double>> &outputClusterMeans, size_t clusterCount, int kMeansIterations)
{
  if (clusterCount == 1) clusterToOne(sigs, clusters, outputClusterMeans);
  else /*if (clusterCount != 100)*/ clusterSignatures(sigs, clusters, outputClusterMeans, clusterCount, kMeansIterations);
  //else tsvqSignatures(sigs, clusters, outputClusterMeans, clusterCount, kMeansIterations);
}

