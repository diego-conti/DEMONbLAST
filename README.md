#DEMONbLAST 0.3

A program to compute **D**iagonal **E**instein **M**etrics **O**n **N**ice **L**ie **A**lgebras of **S**urjective **T**ype.
Based on the following papers:

* D\. Conti and F. A. Rossi, Construction of nice nilpotent Lie groups, [Journal of Algebra, Volume 525, 2019, Pages 311-340, ISSN 0021-8693]( 	
https://doi.org/10.1016/j.jalgebra.2019.01.020).
* D\. Conti and F. A. Rossi, Indefinite Einstein metrics on nice Lie groups, [Forum Mathematicum, vol. 32, no. 6, 2020, pp. 1599-1619]( 	
https://doi.org/10.1515/forum-2020-0049).
* D\. Conti and F. A. Rossi, Ricci-flat and Einstein pseudoriemannian nilmanifolds,  [Complex Manifolds 2019; 6:170-193](	
https://doi.org/10.1515/coma-2019-0010).
* D\. Conti and F. A. Rossi, Nice pseudo-Riemannian nilsolitons. [Journal of Geometry and Physics 173:104433 (2022)](https://doi.org/10.48550/arXiv.2107.07767)

Copyright Diego Conti 2018-2022 diego.conti@unimib.it

###How to build
In order to build *DEMONbLAST*, you need [cmake](https://cmake.org/) and [Wedge 0.4](https://github.com/diego-conti/wedge). Run

	mkdir build
	cd build
	cmake ..
	cmake --build .
	
This will create the executable `demonblast` in the subdirectory build. You should run it from the root directory by

	build/demonblast
This will allow *DEMONbLAST* to employ the structure constants stored in the directory `coefficients`. This is necessary to ensure that the resulting output is consistent through different runs and with the literature quoted above.
	
*DEMONbLAST* has been developed and tested under Ubuntu Linux 21.10 with gcc 11.2.0
	
###How to test

You can run some tests by running

	cd build/test
	ctest .
	




