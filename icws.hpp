#ifndef SKETCH_ICWS
#define SKETCH_ICWS

#include <cmath>
#include <vector>
#include <cstdint>
#include <random>
#include <limits>

#include "hashing.hpp"

using namespace std;

// see Ioffe, Sergey. "Improved consistent sampling, weighted minhash and l1 sketching." Data Mining (ICDM), 2010 IEEE 10th International Conference on. IEEE, 2010.
class ICWS {

	private:
		TabulationHashFunction T1, T2, T3, T4, T5;

	public:
		ICWS(mt19937_64& rng) : T1(rng), T2(rng), T3(rng), T4(rng), T5(rng) {}

		pair<uint64_t, double> operator()(const vector<pair<uint64_t, double>>& x) {
		
			double a_star = numeric_limits<double>::max();
			uint64_t k_star = 0;
			double y_star = 0;

			for(const pair<uint64_t, double>& element : x) {
				uint64_t k = element.first;
				double S_k = element.second;
				double r_u1 = to_unit(T1(k));
				double r_u2 = to_unit(T2(k));
				double c_u1 = to_unit(T3(k));
				double c_u2 = to_unit(T4(k));
				double beta_k = to_unit(T5(k));
				double r_k = -log(r_u1) -log(r_u2); 
				double c_k = -log(c_u1) -log(c_u2); 

				double t_k = floor(log(S_k)/r_k + beta_k);
				double y_k = exp(r_k*(t_k - beta_k));
				double a_k = c_k/(y_k*exp(r_k));
				if(a_k < a_star) {
					a_star = a_k;
					k_star = k;
					y_star = y_k;
				}
			}
			
			return pair<uint64_t, double>(k_star, y_star);
		}
};

class ICWS_t {
	private:
		uint64_t t;
		vector<ICWS> M;
		
	public:
		ICWS_t(mt19937_64& rng, uint64_t t) : t(t) {
			for(uint64_t i = 0; i < t; i++) {
				M.push_back(ICWS(rng));
			}
		}
		
		vector<pair<uint64_t, double>> operator()(const vector<pair<uint64_t, double>>& x) {
			vector<pair<uint64_t, double>> minhashes;
			for(uint64_t i = 0; i < t; i++) {
				minhashes.push_back((M[i])(x));
			}
			return minhashes;
		}
};

// ICWS with tabulated random draws and precomputed logarithms of weights
// Also working with the logarithm og y_k and a_k as suggested in the ICWS paper
// Completely avoids logarithms, exponentials, and divisions. Relies entirely on table lookups and multiplication.
class FastICWS {

	private:
		TabulationHashFunction T1, T2, T3;
		const vector<double>& gamma;
		const vector<double>& gamma_inv;
		const vector<double>& log_gamma;

	public:
		FastICWS(mt19937_64& rng, const vector<double>& gamma, const vector<double>& gamma_inv, const vector<double>& log_gamma) 
			: T1(rng), T2(rng), T3(rng), gamma(gamma), gamma_inv(gamma_inv), log_gamma(log_gamma) {}

		pair<uint64_t, double> operator()(const vector<pair<uint64_t, double>>& log_weight_x) {
		
			const uint64_t MASK16 = 0xFFFFull;
			double log_a_star = numeric_limits<double>::max();
			uint64_t minhash = 0;

			for(const pair<uint64_t, double>& element : log_weight_x) {
				uint64_t k = element.first;
				double log_S_k = element.second;

				uint64_t z = T1(k);
				double r_k = gamma[z & MASK16]; 
				double r_k_inv = gamma_inv[z & MASK16]; 
				double log_c_k = log_gamma[(z >> 16) & MASK16];
				double beta_k = to_unit32(z >> 32);

				double t_k = floor(log_S_k*r_k_inv + beta_k);
				double log_y_k = r_k*(t_k - beta_k);

				double log_a_k = log_c_k - log_y_k - r_k;
				if(log_a_k < log_a_star) {
					log_a_star = log_a_k;
					minhash = T2((uint64_t)t_k) ^ T3(k); // 64-bit minhash
				}
			}
			
			return pair<uint64_t, double>(minhash, log_a_star);
		}
};

class FastICWS_t {
	private:
		uint64_t t;
		vector<FastICWS> M;
		vector<double> gamma;
		vector<double> gamma_inv;
		vector<double> log_gamma;
		
	public:
		FastICWS_t(mt19937_64& rng, uint64_t t) : t(t) {

			// Create discretized version of the X ~ Gamma(2,1) distribution
			// pdf: z*exp(-z)
			// cdf: 1 - exp(-z)*(z + 1)
			// We want to create tables with 2^16 entries
			// The ith entry (starting from 0) will be an interpolation between a value z_{i} with Pr[X <= z_{i}] <= (i + 1)*epsilon
			// and a value z_{i+1} > z_{i} with Pr[X <= z_{i+1}] <= (i+2)*epsilon

			double z = 0.0;
			double z_prev = 0.0;
			double epsilon = 1.0/((1 << 16) + 1);
			double delta = epsilon; // How much to advance z in each step. We can set this to epsilon without skipping steps because pdf <= 1/e.

			for(int i = 0; i < (1 << 16); i++) {
				double target_mass = (i+1)*epsilon;
				while(1 - exp(-z)*(z+1) < target_mass) {
					z += delta;
				}
				// Fill tables
				double v = (z + z_prev)/2;
				z_prev = z;
				gamma.push_back(v);
				gamma_inv.push_back(1/v);
				log_gamma.push_back(log(v));
			}

			for(uint64_t i = 0; i < t; i++) {
				M.push_back(FastICWS(rng, gamma, gamma_inv, log_gamma));
			}
		}
		
		vector<pair<uint64_t, double>> operator()(const vector<pair<uint64_t, double>>& x) {
			vector<pair<uint64_t, double>> log_weight_x;
			for(auto element : x) {
				log_weight_x.push_back({element.first, log(element.second)});
			}
			vector<pair<uint64_t, double>> minhashes;
			for(uint64_t i = 0; i < t; i++) {
				minhashes.push_back((M[i])(log_weight_x));
			}
			return minhashes;
		}
};


#endif
