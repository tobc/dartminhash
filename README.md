# DartMinHash: Fast Sketching for Weighted Sets

This repository contains experiments for comparing the estimation accuracy and running times of the following weighted minwise hashing algorithms:

* DartMinHash https://arxiv.org/abs/2005.11547 
* ICWS and FastICWS https://research.google/pubs/pub36928/
* BagMinHash https://arxiv.org/abs/1802.03914

For BagMinHash and ICWS we use the implementation from here https://github.com/oertl/bagminhash with the relevant files included in the /bagminhash folder.

See the DartMinHash paper https://arxiv.org/abs/2005.11547 for a description of the algorithm and further results of experiments.

## Requirements
The BagMinHash algorithm uses XXHash64 which must be installed:

1. Get xxhash from https://github.com/Cyan4973/xxHash, e.g. using `git clone https://github.com/Cyan4973/xxHash.git`
2. Build using `make lib`
3. Place xxhash.h and libxxhash.a into the directory /bagminhash/xxhash

The code compiles under GCC version 7.5.0 https://gcc.gnu.org/ with relevant commands in the makefile https://www.gnu.org/software/make/.

## Commands

`make run` compiles and executes the main function in main.cpp.

`make test` compiles and run unit tests.

## Experiments
The different experiments are all placed in the main.cpp and write their output to stdout in CSV format. 

1. `time_performance`: Times different algorithms on synthetic data for all combinations of sketch lengths, and L0 and L1 norms chosen. 
2. `time_performance_specific`: Same as above, but only runs on specified tuples of parameters.
3. `measure_similarity`: Returns the estimated Jaccard similarity of different algorithms on synthetic pairs of weighted sets with a specific similarity.

By default `make run` will run `time_performance_specific` on a subset of the settings used in Table 1 in the paper.

In order to pipe the output to the file `data.csv` use command `make run > data.csv`.

## Example output

Notation:

* t denotes the sketch length (usually k in the paper).
* ICWS is a simple and unoptimized version of ICWS using tabulation hashing.
* ICWS_xxhash is the implementation from the BagMinHash repository which uses the ziggurat algorithm for fast sampling: https://en.wikipedia.org/wiki/Ziggurat_algorithm
* FastICWS is our own highly optimized implementation of ICWS that tabulates expensive operations and only computes the logarithms of weights once.
* BagMinHash1 and BagMinHash2: BagMinHash variants described in the BagMinHash paper. BagMinHash2 is essentially always faster and is what we compare against.
* DartMinHash: Optimized implementation following the pseudocode in the paper.

### Performance timings

| id | L0   | log2_L1 | t    | ICWS    | FastICWS | ICWS_xxhash | BagMinHash1 | BagMinHash2 | DartMinHash |
|----|------|---------|------|---------|----------|-------------|-------------|-------------|-------------|
| 0  | 64   | 0.000   | 64   | 0.899   | 0.060    | 0.538       | 2.439       | 0.628       | 0.042       |
| 1  | 1024 | 0.000   | 64   | 11.565  | 0.515    | 9.604       | 4.374       | 1.706       | 0.145       |
| 2  | 64   | 0.000   | 1024 | 19.296  | 2.885    | 8.083       | 48.248      | 13.279      | 0.592       |
| 3  | 1024 | 0.000   | 1024 | 187.661 | 12.643   | 120.135     | 79.775      | 16.586      | 0.824       |
| 4  | 256  | 0.000   | 1    | 0.040   | 0.008    | 0.040       | 0.112       | 0.103       | 0.021       |
| 5  | 256  | 0.000   | 256  | 14.645  | 0.939    | 7.716       | 13.687      | 3.270       | 0.187       |
| 6  | 1024 | 0.000   | 256  | 45.239  | 2.703    | 30.127      | 18.175      | 4.296       | 0.274       |
| 7  | 1024 | 64.000  | 256  | 46.717  | 2.720    | 30.122      | 18.241      | 4.250       | 2.632       |
| 8  | 1024 | -64.000 | 256  | 47.677  | 2.719    | 30.117      | 18.096      | 4.192       | 2.333       |

### Jaccard similarity estimates

| sim_j | t  | ICWS_xxhash | FastICWS | BagMinHash2 | DartMinHash |
|-------|----|-------------|----------|-------------|-------------|
| 0.500 | 1  | 1.000       | 1.000    | 0.000       | 1.000       |
| 0.500 | 2  | 0.500       | 0.500    | 0.000       | 0.500       |
| 0.500 | 3  | 0.333       | 0.333    | 0.000       | 0.333       |
| 0.500 | 4  | 0.500       | 0.250    | 0.750       | 0.750       |
| 0.500 | 5  | 0.000       | 0.400    | 0.600       | 0.200       |
| 0.500 | 6  | 0.667       | 0.500    | 0.500       | 0.000       |
| 0.500 | 7  | 0.571       | 0.714    | 0.429       | 0.429       |
| 0.500 | 8  | 0.250       | 0.375    | 0.625       | 0.500       |
| 0.500 | 9  | 0.889       | 0.222    | 0.556       | 0.444       |
| 0.500 | 10 | 0.600       | 0.400    | 0.700       | 0.400       |

## Tests
We use Catch2 https://github.com/catchorg/Catch2 for unit testing.

To compile and run tests use the command: `make test` 