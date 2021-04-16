#include <cstdio>
#include <iostream>
#include <random>
#include <string>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/SVD>
 
using namespace Eigen;
using namespace std;

class KmerSignatures {
  private:
  unordered_map<string, int> kmerVocab;
  size_t kmerLength;
  size_t sigDim;
  vector<vector<double>> kmerSigs;
  
  vector<vector<double>> getSvd(int dim) const
  {
    default_random_engine rng;
    std::normal_distribution<double> dist(0.0,1.0);
    MatrixXd m(dim, dim);
    for (int y = 0; y < dim; y++) {
      for (int x = 0; x < dim; x++) {
        m(y, x) = dist(rng);
      }
    }
    BDCSVD<MatrixXd> svd;
    svd.compute(m, ComputeFullV);
    auto outputMat = svd.matrixV();
    
    vector<vector<double>> sigs(dim, vector<double>(dim));
    for (int y = 0; y < dim; y++) {
      for (int x = 0; x < dim; x++) {
        sigs[y][x] = outputMat(y, x);
      }
    }
    return sigs;
  }

  string reverseComplement(string input) const
  {
    string output;
    for (size_t i = 0; i < input.size(); i++) {
      char c = input[input.size() - 1 - i];
      switch (c) {
        case 'A': c = 'T'; break;
        case 'C': c = 'G'; break;
        case 'G': c = 'C'; break;
        case 'T': c = 'A'; break;
      }
      output.push_back(c);
    }
    return output;
  }

  string getCanonical(string input) const
  {
    return min(input, reverseComplement(input));
  }

  void getKmerVocab(string prefix, size_t k, set<string> &kmers) const
  {
    static const char nucleotides[] = {'A', 'C', 'G', 'T'};
    for (char n : nucleotides) {
      string kmer = prefix + n;
      if (kmer.size() == k) {
        kmers.insert(getCanonical(kmer));
      } else {
        getKmerVocab(kmer, k, kmers);
      }
    }
  }

  unordered_map<string, int> getKmerVocab(size_t k, size_t &count) const
  {
    unordered_map<string, int> kmerIds;
    set<string> kmers;
    getKmerVocab("", k, kmers);
    int kmerId = 0;
    for (const auto &kmer : kmers) {
      //printf("%s %d\n", kmer.c_str(), kmerId);
      kmerIds[kmer] = kmerId;
      kmerIds[reverseComplement(kmer)] = kmerId;
      kmerId++;
    }
    count = kmerId;
    return kmerIds;
  }
  
  const double *getKmerSignature(const string &kmer) const {
    auto kmerPos = kmerVocab.find(kmer);
    if (kmerPos == kmerVocab.end()) return nullptr;
    int kmerId = kmerPos->second;
    return kmerSigs[kmerId].data();
  }
  
  void normalise(vector<double> &docSig) const
  {
    double total = 0.0;
    for (double v : docSig) {
      total += v * v;
    }
    total = sqrt(total);
    if (total != 0.0) {
      for (double &v : docSig) {
        v /= total;
      }
    }
  }

  public:
  KmerSignatures() {}
  KmerSignatures(size_t k, size_t sigSize) {
    kmerLength = k;
    kmerVocab = getKmerVocab(kmerLength, sigDim);
    if (sigSize == 0) sigSize = sigDim;
    else sigDim = sigSize;
    kmerSigs = getSvd(sigSize);
  }
  
  void getDocumentSignature(const string &seq, vector<double> &docSig) const {
    for (size_t i = 0; i < seq.size() - kmerLength + 1; i++) {
      string kmer = seq.substr(i, kmerLength);
      const double *sig = getKmerSignature(kmer);
      if (sig) {
        //fprintf(stderr, "%s (added)\n", kmer.c_str());
        for (size_t j = 0; j < sigDim; j++) {
          docSig[j] += sig[j];
        }
      } else {
        //fprintf(stderr, "%s (ignored)\n", kmer.c_str());
      }
    }
    normalise(docSig);
  }
  
  size_t getDimensionality() const {
    return sigDim;
  }
  
  size_t getBinSignatureSize() const {
    return (sigDim + 63) / 64;
  }
};
