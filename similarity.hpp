#ifndef SKETCH_SIMILARITY
#define SKETCH_SIMILARITY

#include <vector>
#include <cstdint>
#include <cmath>

using namespace std;

double weight(const vector<pair<uint64_t, double>>&  x) {
    double w = 0;
    for(const pair<uint64_t, double>& v : x) {
        w += v.second;
    }
    return w;
}

double intersection(const vector<pair<uint64_t, double>>& x, const vector<pair<uint64_t, double>>& y) {
	uint64_t i = 0;
	uint64_t j = 0;
	double s = 0;
	while(i < x.size() && j < y.size()) {
		if(x[i].first == y[j].first) {
			s += min(x[i].second, y[j].second);
			i++;
			j++;
		} else if(x[i].first < y[j].first) {
			i++;
		} else {
			j++;
		}
	}
	return s;
}

double jaccard_similarity(const vector<pair<uint64_t, double>>& x, const vector<pair<uint64_t, double>>& y) {
	double s = intersection(x, y);
	double w_x = weight(x);
	double w_y = weight(y);
	return s/(w_x + w_y - s);
}

double l1_similarity(const vector<pair<uint64_t, double>>& x, const vector<pair<uint64_t, double>>& y) {
	double s = intersection(x, y);
	double w_x = weight(x);
	double w_y = weight(y);
	return s/min(w_x, w_y);
}

double hamming_distance(const vector<bool>& x, const vector<bool>& y) {
	double h = 0;
	for(uint32_t i = 0; i < x.size(); i++) {
		if(x[i] != y[i]) {
			h = h + 1;
		}
	}
	return h;
}

double onebit_minhash_jaccard_estimate(const vector<bool>& x, const vector<bool>& y) {
	double h = hamming_distance(x, y);
	double t = x.size();
	return max(0.0, 2*(1 - h/t) - 1);
}

// Similarity conversions
// L1 similarity is the normalized intersection: |x \cap y| / min(|x|, |y|)
// Jaccard similarity is: |x \cap y| / |x \cup y|
double jaccard_similarity_from_l1_similarity(double x_weight, double y_weight, double l1_sim) {
  double i = min(x_weight, y_weight)*l1_sim;
  double u = x_weight + y_weight - i;
  return i/u;
}

double l1_similarity_from_jaccard_similarity(double x_weight, double y_weight, double jaccard_sim) {
  double i = jaccard_sim*(x_weight + y_weight)/(1 + jaccard_sim);
  return i/min(x_weight, y_weight);
}

// Count the number of collisions in two vectors of minhash sketches ((id, rank) pairs)
uint64_t count_collisions(const vector<pair<uint64_t, double>>& x, const vector<pair<uint64_t, double>>& y) {
	uint64_t collisions = 0;
	for(uint64_t i = 0; i < x.size(); i++) {
		if(x[i].first == y[i].first) {
			collisions++;
		}
	}
	return collisions;
}

double jaccard_estimate_from_minhashes(const vector<pair<uint64_t, double>>& x, const vector<pair<uint64_t, double>>& y) {
	return (double)count_collisions(x, y)/x.size();
}

#endif