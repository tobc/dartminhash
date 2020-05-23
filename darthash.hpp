#ifndef SKETCH_DARTHASH
#define SKETCH_DARTHASH

#include <cstdint>
#include <random>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <cmath>
#include "similarity.hpp"
#include "hashing.hpp"

using namespace std;

class DartHash {

	private:
		uint64_t t;
        TabulationHashFunction32 T_nu, T_rho, T_w, T_r;
        TabulationHashFunction T_i, T_p, T_q, F, M;
        vector<double> powers_of_two;
        vector<double> negative_powers_of_two;
        vector<double> poisson_cdf;
        
	public:
		DartHash(mt19937_64& rng, uint64_t t) : t(t), T_nu(rng), T_rho(rng), T_w(rng), T_r(rng), T_i(rng), T_p(rng), T_q(rng), F(rng), M(rng) {
            
            // Tabulate positive and negative powers of two
            double p = 1.0;
            double q = 1.0;
            // Standard double precision ranges from around +-2^1024. We will support slightly less. 
            for(uint64_t i = 0; i < 1000; i++) {
                powers_of_two.push_back(p);
                p = 2.0*p;
                negative_powers_of_two.push_back(q);
                q = 0.5*q;
            }

            // Tabulate the Poisson CDF
            double pdf = exp(-1.0);
            double cdf = pdf;
            for(uint64_t i = 0; i < 100; i++) {
                poisson_cdf.push_back(cdf);
                pdf = pdf/(i + 1);
                cdf += pdf;
            }
        };

        vector<pair<uint64_t, double>> operator()(const vector<pair<uint64_t, double>>& x, double theta = 1.0) {
            vector<pair<uint64_t, double>> darts;
            darts.reserve(2*t);
            double max_rank = theta/weight(x);
            double t_inv = 1.0/t;
            uint32_t RHO = (uint32_t)floor(log2(1.0 + max_rank));

            for(const pair<uint64_t, double>& element : x) {
                uint64_t i = element.first;
                double xi = element.second;
                uint64_t i_hash = T_i(i);
                uint32_t NU = (uint32_t)floor(log2(1.0 + t*xi));
                for(uint32_t nu = 0; nu <= NU; nu++) {
                    uint64_t nu_hash = T_nu(nu);
                    for(uint32_t rho = 0; rho <= RHO; rho++) {
                        uint64_t region_hash = nu_hash ^ T_rho(rho);
                        double two_nu = powers_of_two[nu];
                        double two_rho = powers_of_two[rho];
                        double W = (two_nu - 1)*t_inv;
                        double R = two_rho - 1;
                        double delta_nu = two_nu*t_inv*negative_powers_of_two[rho];
                        double delta_rho = two_rho*negative_powers_of_two[nu];
                        double w0 = W;
                        uint32_t w_max = rho < 32 ? 1ul << rho : 1ul << 31;
                        for(uint32_t w = 0; w < w_max; w++) {
                            if(xi < w0) break;
                            uint64_t w_hash = T_w(w);
                            double r0 = R;
                            uint32_t r_max = nu < 32 ? 1ul << nu : 1ul << 31;
                            for(uint32_t r = 0; r < r_max; r++) {
                                if(max_rank < r0) break;
                                // Get area fingerprint to speed up subsequent hashing 
                                uint64_t area_hash = w_hash ^ T_r(r);
                                uint64_t z = i_hash ^ region_hash ^ area_hash;
                                
                                // Draw from Poisson distribution
                                double p_z = to_unit(T_p(z));
                                uint8_t p = 0;
                                while(p_z > poisson_cdf[p]) {
                                    p++;
                                }
                                
                                uint64_t q = 0;
                                while(q < p) {
                                    // Layer the q-values over z to create a unique key with a strong hash value
                                    uint64_t z_q = z ^ (q << 56) ^ (q << 48) ^ (q << 40) ^ (q << 32) ^ (q << 24) ^ (q << 16) ^ (q << 8) ^ q;
                                    auto uniform_weight_rank = to_units(T_q(z_q));
                                    double weight = w0 + delta_nu*uniform_weight_rank.first;
                                    double rank = r0 + delta_rho*uniform_weight_rank.second;
                                    if(weight < xi && rank < max_rank) {
                                        darts.push_back({F(z_q), rank});
                                    }
                                    q++;
                                }

                                r0 += delta_rho;
                            }
                            w0 += delta_nu;
                        }
                    }
                }
            }
            return darts;
        }


        // Convert the t darts to k minhashes by hashing the darts to k buckets and keeping the minimum from each bucket
        vector<pair<uint64_t, double>> minhash(const vector<pair<uint64_t, double>>& x, uint64_t k) { 
            auto darts = (*this)(x);
            vector<pair<uint64_t, double>> minhashes(k, {0, numeric_limits<double>::max()});
            for(auto& dart : darts) {
                uint64_t j = M(dart.first) % k;
                if(dart.second < minhashes[j].second) {
                    minhashes[j] = dart;
                }
            }
            return minhashes;
        }

        vector<bool> onebit_minhash(const vector<pair<uint64_t, double>>& x, uint64_t k) { 
            vector<bool> sketch(k, false);
            auto minhashes = minhash(x, k);
            for(uint64_t i = 0; i < k; i++) {
                sketch[i] = ((minhashes[i].first & 1ull) == 1ull); // use first bit of MinHash id  
            }
            return sketch;
        }
};

#endif