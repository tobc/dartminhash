#include<vector>
#include<random>
#include <iostream>
#include "catch.hpp"
#include "datagenerator.hpp"
#include "similarity.hpp"
#include "darthash.hpp"
#include "icws.hpp"
#include "dartminhash.hpp"
#include "bagminhash_wrappers.hpp"

using namespace std;

TEST_CASE("Randomly generated weighted sets behave as expected", "[datagenerator]") {
    SECTION("Size and weight") {
        uint64_t seed = 1;
        mt19937_64 rng(seed);
        uint64_t L0 = 128;
        double L1 = 1.0;
        vector<pair<uint64_t, double>> x = generate_weighted_set(L0, L1, rng);
        REQUIRE(x.size() == L0);

        bool sorted = true;
        for(uint32_t i = 1; i < x.size(); i++) {
            if(x[i-1].first > x[i].first) {
                sorted = false;
            }
        }
        REQUIRE(sorted);
        REQUIRE(weight(x) == Approx(L1));
    }
    SECTION("Similar sets") {
        uint64_t seed = 1;
        mt19937_64 rng(seed);
        uint64_t L0 = 128;
        double L1 = 1.0;
        vector<pair<uint64_t, double>> x = generate_weighted_set(L0, L1, rng);
        auto y = generate_similar_weighted_set(x, 0.5, rng);

        // Ensure that the weight is maintained and that we have added an additional element for the excess weight.
        REQUIRE(weight(y) == Approx(L1));
        REQUIRE(y.size() == L0 + 1);

        // Ensure that the generated set matches the desired similarity
        int k = 10;
        for(int i = 0; i <= k; i++) {
            double s = ((double)i/k);
            auto y = generate_similar_weighted_set(x, s, rng);
            REQUIRE(jaccard_similarity(x, y) == Approx(s/(2-s)));
        }
    }
}

TEST_CASE("Similarity measures", "[similarity]") {
    SECTION("weight") {
        vector<pair<uint64_t, double>> x = {
            {1, 1.0},
            {2, 2.0}
        };
        REQUIRE(weight(x) == Approx(3.0));
    }
    SECTION("jaccard_similarity") {
        vector<pair<uint64_t, double>> x = {
            {1, 1.0},
            {2, 2.0}
        };
        vector<pair<uint64_t, double>> y = {
            {1, 1.0},
            {2, 1.0},
            {3, 1.0}
        };
        REQUIRE(intersection(x, y) == Approx(2.0));
        REQUIRE(jaccard_similarity(x, y) == Approx(0.5));
    }
    SECTION("Hamming distance") {
        REQUIRE(hamming_distance({true, true, true}, {false, false, false}) == 3);
        REQUIRE(hamming_distance({true, true, true}, {false, true, false}) == 2);
    }

    SECTION("similarity conversions") {
        vector<pair<uint64_t, double>> x = {
            {1, 1.0},
            {2, 3.0}
        };
        vector<pair<uint64_t, double>> y = {
            {1, 1.0},
            {2, 1.0},
            {3, 1.0}
        };
        // jaccard similarity: 2/5
        // l1 similarity: 2/3
        double jaccard_sim = jaccard_similarity(x, y);
        double l1_sim = l1_similarity(x, y);
        REQUIRE(jaccard_sim == Approx(2.0/5));
        REQUIRE(l1_sim == Approx(2.0/3));
        double x_weight = weight(x);
        double y_weight = weight(y);
        REQUIRE(l1_similarity_from_jaccard_similarity(x_weight, y_weight, jaccard_sim) == Approx(l1_sim));
        REQUIRE(jaccard_similarity_from_l1_similarity(x_weight, y_weight, l1_sim) == Approx(jaccard_sim));
    }
}

TEST_CASE("DartHash", "[darthash]") {
    SECTION("Basic dart properties") {

        uint64_t seed = 1;
        mt19937_64 rng(seed);
        uint64_t t = 256;
	    DartHash D(rng, t);
        uint64_t L0 = 128;
        double L1 = 1.0;
        auto x = generate_weighted_set(L0, L1, rng);
        auto darts = D(x);

        // Number of darts
        // According to http://www.cs.columbia.edu/~ccanonne/files/misc/2017-poissonconcentration.pdf
        // The probability that the number of darts deviates by more than t/2 is at most 2*exp(-t/10)
        REQUIRE(darts.size() > 128);
        REQUIRE(darts.size() < 256 + 128);

        // Dart ranks should be smaller than 1/L1 and fingerprints should be unique
        unordered_set<uint64_t> fingerprints;
        bool too_large = false;
        bool too_small = true;
        for(auto& element : darts) {
            fingerprints.insert(element.first);
            if(element.second > 1/L1) {
                too_large = true;
            }
            
            // Ranks should be uniformly distributed between zero and 1/L1.
            if(element.second > 1/(2*L1)) {
                too_small = false;
            }
        }
        REQUIRE(!too_large);
        REQUIRE(!too_small);
        REQUIRE(fingerprints.size() == darts.size());
    }

    SECTION("Darts to MinHash") {
        // Converting t darts to k minhashes
        // The probability of an "empty" minhash is at most t*exp(-t/k) by a standard union bound over Poisson distributions
        // Verify that when t/k is large that we have no empty minhashes
        // We set the id of empty minhashes to 0
        uint64_t seed = 1;
        mt19937_64 rng(seed);
        uint64_t t = 4096;
        uint64_t k = 128;
	    DartHash D(rng, t);
        uint64_t L0 = 128;
        double L1 = 1.0;
        auto x = generate_weighted_set(L0, L1, rng);
        auto minhashes = D.minhash(x, k);
        bool all_nonempty = true;
        for(auto mh : minhashes) {
            if(mh.first == 0) {
                all_nonempty = false;
            }
        }
        REQUIRE(all_nonempty);

        // When k = t then we expect empty minhashes
        k = t;
        minhashes = D.minhash(x, k);
        all_nonempty = true;
        for(auto mh : minhashes) {
            if(mh.first == 0) {
                all_nonempty = false;
            }
        }
        REQUIRE(!all_nonempty);
    }

    SECTION("MAE of 1-bit minhash sketch stays within Hoeffding bounds") {
        uint64_t seed = 1;
        mt19937_64 rng(seed);
        uint64_t t = 512;
        uint64_t k = 64;
	    DartHash D(rng, t);
        uint64_t L0 = 64;
        double L1 = 1.0;
        double l1_sim = 0.5;
        uint64_t n = 2000;

        // MAE when using minhash to estimate l1 similarity
        double target_l1_sim_mae = 0.1079063;
        double epsilon = 0.05;

        double total_absolute_error = 0.0;
        for(uint64_t i = 0; i < n; i++) {
            auto x = generate_weighted_set(L0, L1, rng);
            auto y = generate_similar_weighted_set(x, l1_sim, rng);
            auto sketch_x = D.onebit_minhash(x, k);
            auto sketch_y = D.onebit_minhash(y, k);
            double jaccard_estimate = onebit_minhash_jaccard_estimate(sketch_x, sketch_y);
            total_absolute_error += abs(l1_similarity_from_jaccard_similarity(weight(x), weight(y), jaccard_estimate) - l1_sim);
        }

        double empirical_mae = total_absolute_error/n;
        REQUIRE(abs(target_l1_sim_mae - empirical_mae) <= epsilon);
    }
}

TEST_CASE("ICWS", "[icws]") {
    SECTION("Weighted samples are valid") {
        uint64_t seed = 1;
        mt19937_64 rng(seed);
        ICWS H(rng);
        int m = 100;
        uint64_t L0 = 64;
        double L1 = 1.0;
        for(int i = 0; i < m; i++) {
            auto x = generate_weighted_set(L0, L1, rng);
            auto z = H(x);
            bool valid_cws = false;
            for(auto element : x) {
                if(element.first == z.first && z.second <= element.second) {
                    valid_cws = true;
                }
            }
            REQUIRE(valid_cws);
        }
    }
}

TEST_CASE("DartMinHash", "[dartminhash]") {
    SECTION("1-bit dartminhash MAE stays within Hoeffding bounds") {

        uint64_t seed = 1;
        mt19937_64 rng(seed);
        uint64_t k = 64;
	    DartMinHash M(rng, k);
        uint64_t L0 = 64;
        double L1 = 1.0;
        double l1_sim = 0.5;
        uint64_t n = 2000;

        // MAE when using minhash to estimate l1 similarity
        double target_l1_sim_mae = 0.1079063;
        double epsilon = 0.05;

        double total_absolute_error = 0.0;
        for(uint64_t i = 0; i < n; i++) {
            auto x = generate_weighted_set(L0, L1, rng);
            auto y = generate_similar_weighted_set(x, l1_sim, rng);
            auto sketch_x = M.onebit_minhash(x);
            auto sketch_y = M.onebit_minhash(y);
            double jaccard_estimate = onebit_minhash_jaccard_estimate(sketch_x, sketch_y);
            total_absolute_error += abs(l1_similarity_from_jaccard_similarity(weight(x), weight(y), jaccard_estimate) - l1_sim);
        }

        double empirical_mae = total_absolute_error/n;
        REQUIRE(abs(target_l1_sim_mae - empirical_mae) <= epsilon);

    }
}

// Test correct estimation of jaccard similarity within Hoeffding bounds
TEST_CASE("Jaccard similarity estimation", "[bagminhash]") {
    uint64_t seed = 1;
    mt19937_64 rng(seed);
    uint64_t L0 = 64;
    double L1 = 1.0;
    double l1_sim = 0.5;
    double target_jaccard_similarity = jaccard_similarity_from_l1_similarity(L1, L1, l1_sim);
    double epsilon = 0.05;
    uint64_t t = 2000;
    auto x = generate_weighted_set(L0, L1, rng);
    auto y = generate_similar_weighted_set(x, l1_sim, rng);
    
    SECTION("BagMinHash1") {
        BagMinHash1 B1(t);
        auto mh_x = B1(x);
        auto mh_y = B1(y);
        double estimated_jaccard_similarity = (double)count_collisions(mh_x, mh_y)/t;
        REQUIRE(abs(target_jaccard_similarity - estimated_jaccard_similarity) <= epsilon);
    }

    SECTION("BagMinHash2") {
        BagMinHash2 B2(t);
        auto mh_x = B2(x);
        auto mh_y = B2(y);
        double estimated_jaccard_similarity = (double)count_collisions(mh_x, mh_y)/t;
        REQUIRE(abs(target_jaccard_similarity - estimated_jaccard_similarity) <= epsilon);
    }

    SECTION("ICWS_xxhash") {
        ICWS_xxhash I(t);
        auto mh_x = I(x);
        auto mh_y = I(y);
        double estimated_jaccard_similarity = (double)count_collisions(mh_x, mh_y)/t;
        REQUIRE(abs(target_jaccard_similarity - estimated_jaccard_similarity) <= epsilon);
    }

    SECTION("FastICWS") {
        FastICWS_t F(rng, t);
        auto mh_x = F(x);
        auto mh_y = F(y);
        double estimated_jaccard_similarity = (double)count_collisions(mh_x, mh_y)/t;
        REQUIRE(abs(target_jaccard_similarity - estimated_jaccard_similarity) <= epsilon);
    }
}


