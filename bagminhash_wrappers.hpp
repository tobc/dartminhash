// Simple wrappers around the bagminhash functions
// Note that we use the floatweightdiscretization as it results in improved performance
#ifndef SKETCH_BAGMINHASH
#define SKETCH_BAGMINHASH

#include <cstdint>
#include <vector>
#include "bagminhash/weighted_minwise_hashing.hpp"

using namespace std;

 vector<pair<uint64_t, double>> weightedhashresult_to_pairs(WeightedHashResult res) {
    vector<pair<uint64_t, double>> output;
    for(uint64_t h : res.hashValues) {
        output.push_back({h, 0.0});
    }
    return output;
}

 vector<tuple<uint64_t, double>> pairs_to_tuples(const vector<pair<uint64_t, double>>& x) {
    vector<tuple<uint64_t, double>> x_tuple;
    for(auto& element : x) {
        x_tuple.push_back(tuple<uint64_t, double>(element.first, element.second));
    }
    return x_tuple;
}

class BagMinHash1 {
    private:
        uint64_t t;
    public:
        BagMinHash1(uint64_t t) : t(t) {};
        vector<pair<uint64_t, double>> operator()(const vector<pair<uint64_t, double>>& x) {
            auto x_tuple = pairs_to_tuples(x);
            WeightedHashResult res = bag_min_hash_1<FloatWeightDiscretization, XXHash64>(x_tuple, t);
            return weightedhashresult_to_pairs(res);

        }
};

class BagMinHash2 {
    private:
        uint64_t t;
    public:
        BagMinHash2(uint64_t t) : t(t) {};
        vector<pair<uint64_t, double>> operator()(const vector<pair<uint64_t, double>>& x) {
            auto x_tuple = pairs_to_tuples(x);
            WeightedHashResult res = bag_min_hash_2<FloatWeightDiscretization, XXHash64>(x_tuple, t);
            return weightedhashresult_to_pairs(res);

        }
};

class ICWS_xxhash {
    private:
        uint64_t t;
    public:
        ICWS_xxhash(uint64_t t) : t(t) {};
        vector<pair<uint64_t, double>> operator()(const vector<pair<uint64_t, double>>& x) {
            auto x_tuple = pairs_to_tuples(x);
            WeightedHashResult res = improved_consistent_weighted_hashing<XXHash64>(x_tuple, t);
            return weightedhashresult_to_pairs(res);
        }
};

#endif