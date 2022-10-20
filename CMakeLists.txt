cmake_minimum_required(VERSION 3.10)
project(demonblast VERSION 0.3)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
add_compile_options(-g -Wctor-dtor-privacy -Wreorder -Wold-style-cast -Wsign-promo -Wchar-subscripts -Winit-self -Wmissing-braces -Wparentheses -Wreturn-type -Wswitch -Wtrigraphs -Wextra -Wno-sign-compare -Wno-narrowing -Wno-attributes)


set(SOURCES_NO_MAIN src/partitions.cpp src/tree.cpp src/labeled_tree.cpp src/weightbasis.cpp src/niceliegroup.cpp src/liegroupsfromdiagram.cpp src/gauss.cpp src/log.cpp src/niceeinsteinliegroup.cpp src/ricci.cpp src/filter.cpp src/permutations.cpp src/weightmatrix.cpp src/implicitmetric.cpp src/antidiagonal.cpp src/adinvariantobstruction.cpp src/parsetree.cpp src/automorphisms.cpp src/partitionprocessor.cpp src/diagramprocessor.cpp src/nicediagramsinpartition.cpp)

set(SOURCES src/nice.cpp ${SOURCES_NO_MAIN})

set (INCLUDES src/arrow.h src/labeled_tree.h src/partitions.h src/liegroupsfromdiagram.h src/ permutations.h src/diagramprocessor.h src/linearinequalities.h src/ricci.h src/double_arrows_tree.h src/linearsolve.h src/taskrunner.h src/filter.h src/log.h src/tree.h src/gauss.h src/niceeinsteinliegroup.h src/weightbasis.h src/horizontal.h src/niceliegroup.h src/weightmatrix.h src/ xginac.h src/tree.hpp matrixbuilder.h src/options.h src/implicitmetric.h src/antidiagonal.h src/nicediagramsinpartition.h src/adinvariantobstruction.h src/includes.h src/diagramanalyzer.h src/parsetree.h src/automorphisms.h src/components.h src/coefficientconfiguration.h src/expressionparser.h src/partitionprocessor.h src/coefficientconfiguration.h)

link_libraries(ginac wedge cocoa gmp cln boost_program_options)
link_directories (/home/diego/wedge/lib)
include_directories(/home/diego/wedge/include)

add_executable(demonblast ${SOURCES})
add_executable(debug EXCLUDE_FROM_ALL ${SOURCES} )
target_compile_options(demonblast PUBLIC -O3)
target_compile_options(debug PUBLIC -O0)

#add targets in subdirectory programs, if it exists. This can be used to build other programs that use the code of demonblast
if (EXISTS ${CMAKE_SOURCE_DIR}/programs)	
	add_subdirectory(programs EXCLUDE_FROM_ALL)
endif()
add_subdirectory(test)



