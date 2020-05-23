#ifndef SKETCH_HASHING
#define SKETCH_HASHING

#include <cstdint>
#include <random>

using namespace std;

 double to_unit(uint64_t x) {
    return (double)x/0xFFFFFFFFFFFFFFFFull;
}

 double to_unit32(uint32_t x) {
	return (double)x/0xFFFFFFFFul;
}

// Convert a 64-bit uint to two doubles by splitting it in two and normalizing each 32-bit part
 pair<double, double> to_units(uint64_t x) {
    return { 
        ((double)(x >> 32))/0xFFFFFFFFul,
        (double)(x & 0xFFFFFFFFull)/0xFFFFFFFFul
    };
}

class TabulationHashFunction {

	private:
		const uint64_t mask = 0xFF; 
		uint64_t T1[256];
		uint64_t T2[256];
		uint64_t T3[256];
		uint64_t T4[256];
        uint64_t T5[256];
		uint64_t T6[256];
		uint64_t T7[256];
		uint64_t T8[256];

	public:
		TabulationHashFunction(mt19937_64& rng) {
			for(int i = 0; i < 256; i++) {
				T1[i] = rng();
				T2[i] = rng();
				T3[i] = rng();
				T4[i] = rng();
                T5[i] = rng();
				T6[i] = rng();
				T7[i] = rng();
				T8[i] = rng();
			}
		}

		uint64_t operator()(uint64_t x) {
			uint64_t hashvalue = T1[x & mask] ^ T2[(x >> 8) & mask] ^ T3[(x >> 16) & mask] ^ T4[(x >> 24) & mask] ^ 
                                 T5[(x >> 32) & mask] ^ T6[(x >> 40) & mask] ^ T7[(x >> 48) & mask] ^ T8[(x >> 56) & mask];
			return hashvalue;
		}
};

class TabulationHashFunction8 {

	private:
		uint64_t T1[256];

	public:
		TabulationHashFunction8(mt19937_64& rng) {
			for(int i = 0; i < 256; i++) {
				T1[i] = rng();
			}
		}

		uint64_t operator()(uint8_t x) {
			return T1[x];
		}
};

class TabulationHashFunction32 {

	private:
		const uint32_t mask = 0xFF; 
		uint64_t T1[256];
		uint64_t T2[256];
		uint64_t T3[256];
		uint64_t T4[256];

	public:
		TabulationHashFunction32(mt19937_64& rng) {
			for(int i = 0; i < 256; i++) {
				T1[i] = rng();
				T2[i] = rng();
				T3[i] = rng();
				T4[i] = rng();
			}
		}

		uint64_t operator()(uint32_t x) {
			uint64_t hashvalue = T1[x & mask] ^ T2[(x >> 8) & mask] ^ T3[(x >> 16) & mask] ^ T4[(x >> 24) & mask];
			return hashvalue;
		}
};

#endif