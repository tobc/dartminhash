# DartMinHash: Fast Sketching for Weighted Sets

This repository contains experiments for comparing the estimation accuracy and running times of the following weighted minwise hashing algorithms:

* DartMinHash 
* ICWS and FastICWS https://research.google/pubs/pub36928/
* BagMinHash https://arxiv.org/abs/1802.03914

For BagMinHash and ICWS we use the implementation from here https://github.com/oertl/bagminhash with the relevant files included in the /bagminhash folder.

## Requirements
The BagMinHash algorithm uses XXHash64 which must be installed:

1. Get xxhash from https://github.com/Cyan4973/xxHash, e.g. using `git clone https://github.com/Cyan4973/xxHash.git`
2. Build using `make lib`
3. Place xxhash.h and libxxhash.a into the directory /bagminhash/xxhash

The code compiles under GCC version 7.5.0 https://gcc.gnu.org/ with relevant commands in the makefile https://www.gnu.org/software/make/.

## Commands

`make run` - Compiles and executes the main function in main.cpp.
`make test` - Compiles and run unit tests.

## Experiments
The different experiments are all placed in the main.cpp and write their output to stdout in CSV format. 

1. `time_performance`: Times different algorithms on synthetic data for all combinations of sketch lengths, and L0 and L1 norms chosen. 
2. `time_performance_specific`: Same as above, but only runs on specified tuples of parameters.
3. `measure_similarity`: Returns the estimated Jaccard similarity of different algorithms on synthetic pairs of weighted sets with a specific similarity.

By default `make run` will run `time_performance_specific` on a subset of the settings used in Table 1 in the paper.
In order to pipe the output to the file `data.csv` use command `make run > data.csv`.

## Example output




## Tests
We use Catch2 https://github.com/catchorg/Catch2 for unit testing.

To compile and run test use the command: `make test`



 