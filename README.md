#DEMONbLAST 0.3

A program to compute **D**iagonal **E**instein **M**etrics **O**n **N**ice **L**ie **A**lgebras of **S**urjective **T**ype.
Based on the following papers:

* \[CR19] D. Conti and F. A. Rossi, Construction of nice nilpotent Lie groups, [Journal of Algebra, Volume 525, 2019, Pages 311-340, ISSN 0021-8693](https://doi.org/10.1016/j.jalgebra.2019.01.020).
* \[CR19b] D. Conti and F. A. Rossi, Ricci-flat and Einstein pseudoriemannian nilmanifolds,  [Complex Manifolds 2019; 6:170-193](
https://doi.org/10.1515/coma-2019-0010).
* \[CR20] D. Conti and F. A. Rossi, Indefinite Einstein metrics on nice Lie groups, [Forum Mathematicum, vol. 32, no. 6, 2020, pp. 1599-1619](https://doi.org/10.1515/forum-2020-0049).
* \[CR22] D. Conti and F. A. Rossi, Nice pseudo-Riemannian nilsolitons. [Journal of Geometry and Physics 173:104433 (2022)](https://doi.org/10.48550/arXiv.2107.07767)

Copyright Diego Conti 2018-2022 diego.conti@unimib.it. Latest version available at [https://github.com/diego-conti/DEMONbLAST](https://github.com/diego-conti/DEMONbLAST).

##How to build
In order to build `demonblast`, you need [cmake](https://cmake.org/) and [Wedge 0.4](https://github.com/diego-conti/wedge). Run

	mkdir build
	cd build
	cmake ..
	cmake --build .
	
This will create the executable `demonblast` in the subdirectory build.
	
`demonblast` has been developed and tested under Ubuntu Linux 21.10 with gcc 11.2.0
	
##How to test

You can run some tests by running

	cd build/test
	ctest .
	
(or `ctest --test-dir build/test` with newer versions of `cmake`)

##Running `demonblast`
You should run `demonblast` from the project root directory, e.g.

	build/demonblast --help

This command will print a list of command-line options, some of which are self-explanatory.

#### Diagram selection
There are three ways to run `demonblast`.

- -`-all-partitions n` will process all nice diagrams in dimension n. For example,

	build/demonblast --all-partitions 6

will classify all nice diagrams of dimension 6 and write the output in the directory `output/6`.

- `--partition n1 n2 ... nk` will process the partition *(n1 n2 ... nk)*. For example,
	
	build/demonblast --partition 4 2

will classify all nice diagrams with partition (4,2) and write the output to the file `output/6/part4_2.dot`

- `--digraph diag` will process the single nice diagram indicate by diag. For example,

	build/demonblast --digraph "6:4->[6]1,5->[6]2"

will process the diagram with six nodes and four arrows corresponding to the Lie algebra (46,56,0,0,0,0). Notice that the input format for this option is the same as the format used in the output (indeed the file `output/6/part4_2.dot` contains this string). Notice that the quotation marks are necessary to avoid the character > from being interpreted from the shell as redirection.

Notice that diagrams are cached in the directory `diagrams`. All paths used by `demonblast` are relative to the directory from which it is invoked, which is why invocation from the project root directory is recommended.

This will allow `demonblast` to employ the structure constants stored in the directory `coefficients`. This is necessary to ensure that the resulting output is consistent through different runs and with the literature quoted above.

#### Lie algebra mode

By default, `demonblast` runs in Lie algebra mode. This means that it attempts to classify all Lie algebras for every diagram. This is only possible when the quadratic equations corresponding to the Jacobi identity reduce to linear equations (which is very often the case in low dimensions, since most structure constants can be assumed to be ±1 by exploiting the action of diagonal matrices, see [CR2019]).

Since GiNaC does not produce canonical output, the resulting structure constants may differ in subsequent runs. Additionally, there are cases in which `demonblast` is not able to solve all equations, and they must be handled manually.

To address this, `demonblast` allows caching the structure constants by means of the option `--coefficients store`. This has the effect of storing the structure constants in files contained in the directory `coefficients`. One can then instruct `demonblast` to retrieve structure constants from these files by means of the option `--coefficients load`. 

Notice that the coefficient files can be edited manually, which allows one to plugin explicit solutions for the quadratic equations, or simply to rename the parameters.

`demonblast` is currently shipped with coefficient files up to dimension 9, consistently with the tables of [CR19].

#### Diagram mode

Invoking `demonblast` with the option `--mode diagram` has the effect of classifying nice diagrams without regard for the underlying Lie algebra structures.

#### Table mode

Invoking `demonblast` with the option `--mode table` has the effect of printing out a table of Lie algebras corresponding to the indicated partitions. 

**BUG**: this is printed to a file in --partition mode.

### Metrics

The option `--diagonal-ricci-flat-metrics` instructs `demonblast` to print out the signatures  of diagonal Ricci-flat metrics as in [CR19b].

The option `--diagonal-nilsoliton-metrics` instructs `demonblast` to print out the signatures of diagonal nilsoliton metrics as in [CR22]. In this context, this means that the Ricci operator takes the form Ric=lambda Id + D, with lambda a nonzero constant and D a derivation. In particular, one recovers diagonal Einstein metrics of nonzero scalar curvature in the case where D is zero. Notice that a nilsoliton metric as above is Einstein if and only if all derivations have zero trace, i.e. the Nikolayevsky derivation is zero. Thus, the diagonal Einstein metrics of [CR20] can be computed by invoking

	build/demonblast  --all-partitions 8 --mode table --diagonal-nilsoliton-metrics --only-traceless-derivations

### Restrictions

The following options have the effect of restricting the diagrams that are considered

- `--only-traceless-derivations` only include diagrams for which the Nikolayevsky derivation is zero
- `--simple-Nikolayevsky true` only include diagrams for which the Nikolayevsky derivation is simple (i.e. has distinct eigenvalues)
- `--simple-Nikolayevsky false` only include diagrams for which the Nikolayevsky derivation is simple (i.e. has distinct eigenvalues)
- `--kernel-root-matrix-dimension arg` only include diagrams for which the dimension of the kernel of the root matrix satisfies the indicated equality or inequality. For instance, 
	`--kernel-root-matrix-dimension ">0"` selects diagrams for which the kernel has positive dimension, and
	`--kernel-root-matrix-dimension "=0"` is equivalent to `--only-traceless-derivations`
- `--cokernel-root-matrix-dimension arg` only include diagrams for which the dimension of the cokernel of the root matrix satisfies the indicated equality or inequality.
- `--irreducible` only include connected nice diagrams, i.e. irreducible nice Lie algebras.
- `--only-with-metric` if either `--diagonal-ricci-flat-metrics` or `--diagonal-nilsoliton-metrics` is specified, omit nice diagrams where the indicated type of metric does not appear because of the linear obstructions.

### Other output

The following options control output:

- `--derivations` 	include Lie algebra derivations in output.
- `--lcs-and-ucs` in table mode, write the dimensions of the lower and upper central series.
- `--invert` inverts the nodes, so that e.g. the Heisenberg Lie algebra appears as (0,0,12) instead of (23,0,0).
- `--matrix-data` prints out information associated to the root matrix.
	