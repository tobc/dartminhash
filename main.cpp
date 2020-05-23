#include <cstdint>
#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>

#include "timer.hpp"
#include "icws.hpp"
#include "darthash.hpp"
#include "dartminhash.hpp"
#include "datagenerator.hpp"
#include "bagminhash_wrappers.hpp"

using namespace std;

template <typename T>
double time_algorithm(const vector<vector<pair<uint64_t, double>>>& data, T& hasher, uint64_t& id_xor, double& rank_sum) {
	Timer timer;
	for(auto& x : data) {
		timer.start();
		auto minhashes = hasher(x);
		timer.stop();

		// Do something with the minhashes to ensure that the compiler doesn't optimize them away
		for(auto mh : minhashes) {
			id_xor = id_xor ^ mh.first;
			rank_sum += mh.second;
		}
	}

	return timer.elapsed_ms()/data.size();
}

void time_performance() {
	uint64_t seed = 1;
	mt19937_64 rng(seed);

	// Settings for L0 plot
	// vector<uint64_t> L0_values;
	// for(uint64_t i = 0; i < 17; i++) {
	// 	L0_values.push_back(1ull << i);
	// }
	// vector<double> L1_values = {1.0};
	// vector<uint64_t> t_values = {256};
	// uint64_t m = 100;

	// Settings for L1 plot
	// vector<uint64_t> L0_values = {256};
	// vector<uint64_t> t_values = {256};
	// vector<double> L1_values;
	// for(int i = -128; i <= 128; i = i + 16) {
	// 	L1_values.push_back(pow(2.0, (double)i));
	// }
	// uint64_t m = 100;

	// Settings for k plot
	// vector<uint64_t> t_values;
	// for(uint64_t i = 0; i < 17; i++) {
	// 	t_values.push_back(1ull << i);
	// }
	// vector<double> L1_values = {1.0};
	// vector<uint64_t> L0_values = {256};
	// uint64_t m = 100;

	// Exploration
	// vector<uint64_t> L0_values = {1, 4, 16, 64, 256, 1024, 4096, 16384};
	// vector<double> L1_values = {pow(2.0, -32), pow(2.0, -16), pow(2.0, -8), pow(2.0, 0), pow(2.0, 8), pow(2.0, 16), pow(2.0, 32)};
	// vector<uint64_t> t_values = {1, 4, 16, 64, 256, 1024, 4096, 16384};
	// uint64_t m = 10;

	// Settings for heatmap 
	vector<uint64_t> L0_values;
	vector<uint64_t> t_values;
	for(uint64_t i = 0; i < 15; i++) {
		L0_values.push_back(1ull << i);
		t_values.push_back(1ull << i);
	}
	vector<double> L1_values = {1.0};
	uint64_t m = 100;


	uint64_t id_xor = 0;
	double rank_sum = 0;
	uint64_t experiment_counter = 0;

	cout << setprecision(5) << fixed;
	cout << "id, L0, log2_L1, t, FastICWS, BagMinHash2, DartMinHash" << endl;

	for(uint64_t L0 : L0_values) {
		for(double L1 : L1_values) {
			for(uint64_t t : t_values) {

				cout << experiment_counter << ", " << L0 <<  ", " << log2(L1) <<  ", "  << t << ", ";

				// Generate data
				vector<vector<pair<uint64_t, double>>> data;
				for(uint64_t i = 0; i < m; i++) {
					data.push_back(generate_weighted_set(L0, L1, rng));
				}

				// Algorithms
				// ICWS_t I(rng, t);
				// cout << time_algorithm<ICWS_t>(data, I, id_xor, rank_sum) << ", ";

				// ICWS_xxhash I_xx(t);
				// cout << time_algorithm<ICWS_xxhash>(data, I_xx, id_xor, rank_sum) << ", ";

				FastICWS_t F(rng, t);
				cout << time_algorithm<FastICWS_t>(data, F, id_xor, rank_sum) << ", ";

				// BagMinHash1 B1(t);
				// cout << time_algorithm<BagMinHash1>(dat a, B1, id_xor, rank_sum) << ", ";
				
				BagMinHash2 B2(t);
				cout << time_algorithm<BagMinHash2>(data, B2, id_xor, rank_sum) << ", ";

				DartMinHash M(rng, t);
				cout << time_algorithm<DartMinHash>(data, M, id_xor, rank_sum) << endl;

				experiment_counter++;
			}
		}
	}

	cout << "rank sum: " << rank_sum << ", id XOR: " << id_xor << endl;
}

struct experiment_settings {
	uint64_t L0;
	double L1;
	uint64_t t;
};

void time_performance_specific() {
	uint64_t seed = 1;
	mt19937_64 rng(seed);

	vector<experiment_settings> settings {
		// L0, L1, t
		// Varying L0
		{64, pow(2.0, 0.0), 64},  
		{1024, pow(2.0, 0.0), 64}, 
		// {16384, pow(2.0, 0.0), 64}, 

		{64, pow(2.0, 0.0), 1024},  
		{1024, pow(2.0, 0.0), 1024}, 
		// {16384, pow(2.0, 0.0), 1024}, 

		// // Varying t
		{256, pow(2.0, 0.0), 1}, 
		{256, pow(2.0, 0.0), 256}, 
		// {256, pow(2.0, 0.0), 4096},

		// {4096, pow(2.0, 0.0), 1}, 
		// {4096, pow(2.0, 0.0), 256}, 
		// {4096, pow(2.0, 0.0), 4096},
		
		// // Varying L1
		{1024, pow(2.0, 0.0), 256}, 
		{1024, pow(2.0, 64.0), 256}, 
		{1024, pow(2.0, -64.0), 256}, 
		// {1024, pow(2.0, 512.0), 256}, 
		// {1024, pow(2.0, -512.0), 256}, 
	};
	uint64_t m = 100;

	uint64_t id_xor = 0;
	double rank_sum = 0;
	uint64_t experiment_counter = 0;

	cout << setprecision(3) << fixed;
	cout << "id, L0, log2_L1, t, ICWS, FastICWS, ICWS_xxhash, BagMinHash1, BagMinHash2, DartMinHash" << endl;

	for(experiment_settings s : settings) {

		cout << experiment_counter << ", " << s.L0 <<  ", " << log2(s.L1) <<  ", "  << s.t << ", ";

		// Generate data
		vector<vector<pair<uint64_t, double>>> data;
		for(uint64_t i = 0; i < m; i++) {
			data.push_back(generate_weighted_set(s.L0, s.L1, rng));
		}
		
		// Algorithms
		ICWS_t I(rng, s.t);
		cout << time_algorithm<ICWS_t>(data, I, id_xor, rank_sum) << ", ";

		FastICWS_t F(rng, s.t);
		cout << time_algorithm<FastICWS_t>(data, F, id_xor, rank_sum) << ", ";

		ICWS_xxhash I_xx(s.t);
		cout << time_algorithm<ICWS_xxhash>(data, I_xx, id_xor, rank_sum) << ", ";

		BagMinHash1 B1(s.t);
		cout << time_algorithm<BagMinHash1>(data, B1, id_xor, rank_sum) << ", ";
		
		BagMinHash2 B2(s.t);
		cout << time_algorithm<BagMinHash2>(data, B2, id_xor, rank_sum) << ", ";

		DartMinHash M(rng, s.t);
		cout << time_algorithm<DartMinHash>(data, M, id_xor, rank_sum) << endl;

		experiment_counter++;
	}

	cout << "rank sum: " << rank_sum << ", id XOR: " << id_xor << endl;
}

// Measure how the estimated jaccard similarity changes as the number of minhashes increases
void measure_similarity() {
	uint64_t seed = 1;
	mt19937_64 rng(seed);

	// Experiment in paper
	// uint64_t L0 = 256;
	// double L1 = 1.0;
	// uint64_t t_min = 1;
	// uint64_t t_max = 100;
	// vector<double> jaccard_similarity_values = {0.25, 0.5, 0.75};

	uint64_t L0 = 256;
	double L1 = 1.0;
	uint64_t t_min = 1;
	uint64_t t_max = 10;
	vector<double> jaccard_similarity_values = {0.5};


	vector<uint64_t> t_values;
	for(uint64_t i = t_min; i < t_max + 1; i++) {
		t_values.push_back(i);
	}

	cout << setprecision(3) << fixed;

	cout << "sim_j, t, ICWS_xxhash, FastICWS, BagMinHash2, DartMinHash" << endl;
	for(double jaccard_similarity : jaccard_similarity_values) {
		for(uint64_t t : t_values) {
			cout << jaccard_similarity << ", " << t << ", ";

			// Generate a pair of similar points	
			double l1_sim = l1_similarity_from_jaccard_similarity(L1, L1, jaccard_similarity);
			auto x = generate_weighted_set(L0, L1, rng);
			auto y = generate_similar_weighted_set(x, l1_sim, rng);

			// // Algorithms
			ICWS_xxhash I_xx(t);
			cout << jaccard_estimate_from_minhashes(I_xx(x), I_xx(y)) << ", ";

			FastICWS_t F(rng, t);
			cout << jaccard_estimate_from_minhashes(F(x), F(y)) << ", ";

			BagMinHash2 B2(t);
			cout << jaccard_estimate_from_minhashes(B2(x), B2(y)) << ", ";

			DartMinHash D(rng, t);
			cout << jaccard_estimate_from_minhashes(D(x), D(y)) << endl;

		}
	}
}

int main() {
	// time_performance();
	time_performance_specific();
	// measure_similarity();
}