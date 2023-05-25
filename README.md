# What is gridCHP++?

GridCHP++ is a stabilizer simulator inspired by Scott Aaronson's CHP,
specialized to simulating random Clifford measurements on the grid
Given an sqrt(n) x sqrt(n) cluster state, it simulates X/Y/Z measurements on all qubits

It implements three variants (with increasing efficiency):
1. Naive algorithm: Creates the entire state in memory, and measures
2. Row-by-row algorithm: Only keeps two rows of the state in memory at a time
3. Recursive algorithm: Only keeps certain recursive subblocks in memory

The recursive algorithm runs in time O(n^(3/2)) and space O(n^(1/2)).
See the details in [Fast simulation of planar Clifford circuits](https://arxiv.org/abs/2009.03218) by David Gosset, Daniel Grier, Alex Kerzner, Luke Schaeffer

# Using gridCHP++

Perhaps unsurprisingly, gridCHP++ is written in C++. Compilation is achieved using the makefile.
By default make will compile using gcc. Use 'make clang' to compile with clang.
(this is important because the compiler flags differ in the two cases)

usage:
    ./gridCHP++ [options] <grid length>
where the options are:
    -b <filename>   Set the measurement bases
    -r <filename>   Set the measurement outcomes
    -c              Print circuit instructions (mimics those used in [Stim](https://github.com/quantumlib/Stim))
    -l              Print *only* the running time
    -t              Display the running time
    -s              Hide (silence) measurement outcomes
    -h              Use row-by-row ("holographic") algorithm
    -n              Use naive brute force algorithm
    -o              Use naive brute force algorithm with random ordering of measurements
    -v              Brute force circuit for verification (only used for debugging)
    -d              Simulation without destabilizers

By default, gridCHP++ uses the fast recursive algorithm.

To do timing tests, run
./stats.sh [step size] <grid length>
which will output a stats.csv file containing the runtime for each grid
starting from 2, incrementing by the step size, and stopping at the grid length
