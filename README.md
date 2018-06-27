# MPI_Verification
Verification of MPI programs using Static Analysis. Two tools are used here which together help to verify deadlock-freedom for MPI programs with point-to-point communication calls. These tools are: parallelCFG and Split.

# parallelCFG
## Installation
1. Setup ROSE Development Branch on your system. Refer to the links http://rosecompiler.org/ROSE_HTML_Reference/installation.html and https://en.wikibooks.org/wiki/ROSE_Compiler_Framework/Installation , whichever looks convenient.
2. Go to `.../rose-develop/projects/symbolicAnalysisFramework/src/` and replace the parallelCFG implementation with the one given in this repository.
3. Now go to `.../build-rose/projects/symbolicAnalysisFramework/src/parallelCFG` , delete everything except Makefile and run the `$ make` command
4. Use the command `$ source set.rose` to add the relevant Rose directories to the PATH variable.

## Usage
The `parallelCFG/tests` folder contains the benchmark programs (with pragmas added). Use the command `$ pcfgiterator <name_of_file>` to run pCFG analysis on the program. The output results are added in a file `matches.txt`

# Split
The SPLIT tool has been added in the repository, so you can directly use it.

## Usage
1. Go to `SPLIT/SPLIT tool` and run the command `$ java -jar Split.jar`
2. GUI interface opens. Open any `.spic` file and run the Safety Verifier to check whether the deadlock-freedom condition is an invariant or not.
