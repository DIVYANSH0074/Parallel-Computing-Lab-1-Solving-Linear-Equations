# Parallel Computing Lab 1: Solving Linear Equations
By Jason Yao, [github](https://github.com/JasonYao/Parallel-Computing-Lab-1-Solving-Linear-Equations)

## FOR THE GRADER: 
This README utilises github's markdown, and is a much easier read on the github website listed above.

Link: https://github.com/JasonYao/Parallel-Computing-Lab-1-Solving-Linear-Equations

For my conclusion, please see [here](CONCLUSION.md)

## Description
This program is designed to solve a group of linear equations in an algorithmic fashion. 
Both parallel and sequential versions of this code are provided in the source, and can easily be toggled.

## Compilation & Running
### To compile the code and test against an input automatically
To change the input, simply change `testNumber="01"` to the test case you'd like: small.txt == 01, medium.txt == 02, large.txt == 03, huge.txt == 04
```sh
./compileAndTest
```

### To compile and run the code manually
```sh
mpicc -g -Wall -o <output_file_name> <source_file_name>
mpirun -n <number_of_processes> ./<output_file_name> <inputfile.txt>
```

e.g.
```sh
mpicc -g -Wall -o gs gs.c
mpirun -n 4 ./gs testing/input/input-01
```

Where:
- `mpicc` is a wrapper scipt to include mpi libraries upon compilation

- `-g` is a flag to produce debugging information

- `-Wall` is a flag to turn on all warnings

- `-o` is to signify the output binary

## Debugging version
If you'd like to run the debugging version of the code irrespective of the other flags, please edit the source file [gs.c](gs.c) and change line **12** from

```sh
bool IS_DEBUG_MODE = false;
```
to
```sh
bool IS_DEBUG_MODE = true;
```
then compile again before running

## Sequential version
If you'd like to run the sequential version of the code instead, please edit the source file [gs.c](gs.c) and change line **13** from

```sh
bool IS_SEQUENTIAL_MODE = false;
```
to
```sh
bool IS_SEQUENTIAL_MODE = true;
```
then compile again before running

## Timed version
If you'd like to run the timed version of this code, irrespective of the other flags, please edit the source file [gs.c](gs.c) and change line **14** from

```sh
bool IS_TIMED_MODE = false;
```
to
```sh
bool IS_TIMED_MODE = true;
```
then compile again before running

## License
This repo is licensed under the terms of the GNU GPL v3, a copy of which may be found [here](LICENSE).
