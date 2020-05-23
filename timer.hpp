
#ifndef SKETCH_TIMER
#define SKETCH_TIMER

#include <iostream>
#include <chrono>
#include <iomanip>

using namespace std;

struct Timer {
	string name;
	chrono::nanoseconds elapsed;
	chrono::high_resolution_clock::time_point start_time; 

	Timer() : Timer("Timer") {}
	Timer(string name) : name(name) {
		elapsed = chrono::nanoseconds(0);
	}

	void start() {
		start_time = std::chrono::high_resolution_clock::now();
	}

	void stop() {
		elapsed += std::chrono::high_resolution_clock::now() - start_time;
	}

	void reset() {
		elapsed = chrono::nanoseconds(0);
	}

	double elapsed_s() {
		return elapsed_ms()/1000;
	}

	double elapsed_ms() {
		return (double)elapsed.count()/1000000;
	}

	double elapsed_ns() {
		return elapsed.count();
	}

	void print_ms() {
		cout << fixed << setprecision(2);
		cout << name << ": " << elapsed_ms() << " (ms)" << endl;
	}

	void print_s() {
		cout << fixed << setprecision(2);
		cout << name << ": " << elapsed_s() << " (s)" << endl;
	}
};

#endif
