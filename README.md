# Parallel Computing Lab 1: Solving Linear Equations
By Jason Yao, [github](https://www.github.com/JasonYao/SolvingLinearEquations)

## Note for the grader's convenience (and removal of eyestrain)
This README utilises github's markdown, and as such is much better read on the github website listed above.

## Description
This program is designed to solve a group of linear equations in an algorithmic fashion. 
Both parallel and sequential versions of this code are provided in the source, and can easily be toggled.

## Compilation & Running
### To compile the code and test against an input automatically
(To change the input, simply change `testNumber="01"` to the test case you'd like: small == 01, medium == 02, large == 03, huge == 04) 
```sh
./compileAndTest
```

### To compile and run the code manually
```sh
mpicc -g -Wall -o <output_file_name> <source_file_name>
mpirun -n <number_of_processes> ./gs <inputfile.txt>
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

- `gs` is the output binary's filename

- `gs.c` is the source file

## Sequential version
If you'd like to run the sequential version of the code instead, please set

```sh
bool IS_SEQUENTIAL_MODE = false;
```
to
```sh
bool IS_SEQUENTIAL_MODE = true;
```
then compile again before running

## License
This repo is licensed under the terms of the GNU GPL v3 license, a copy of which may be found [here](LICENSE).
