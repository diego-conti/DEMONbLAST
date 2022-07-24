cmake_minimum_required(VERSION 3.10)
project(demonblast VERSION 0.2)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
add_compile_options(-g -Wctor-dtor-privacy -Wreorder -Wold-style-cast -Wsign-promo -Wchar-subscripts -Winit-self -Wmissing-braces -Wparentheses -Wreturn-type -Wswitch -Wtrigraphs -Wextra -Wno-sign-compare -Wno-narrowing -Wno-attributes)


set(SOURCES_NO_MAIN partitions.cpp tree.cpp labeled_tree.cpp weightbasis.cpp niceliegroup.cpp liegroupsfromdiagram.cpp gauss.cpp log.cpp niceeinsteinliegroup.cpp ricci.cpp filter.cpp permutations.cpp weightmatrix.cpp implicitmetric.cpp antidiagonal.cpp adinvariantobstruction.cpp parsetree.cpp automorphisms.cpp partitionprocessor.cpp diagramprocessor.cpp nicediagramsinpartition.cpp)

set(SOURCES nice.cpp ${SOURCES_NO_MAIN})

set (INCLUDES arrow.h labeled_tree.h partitions.h liegroupsfromdiagram.h  permutations.h diagramprocessor.h linearinequalities.h ricci.h double_arrows_tree.h linearsolve.h taskrunner.h filter.h log.h tree.h gauss.h niceeinsteinliegroup.h weightbasis.h horizontal.h niceliegroup.h weightmatrix.h  xginac.h tree.hpp matrixbuilder.h options.h implicitmetric.h antidiagonal.h nicediagramsinpartition.h adinvariantobstruction.h includes.h diagramanalyzer.h parsetree.h automorphisms.h components.h coefficientconfiguration.h expressionparser.h partitionprocessor.h coefficientconfiguration.h)

link_libraries(ginac wedge cocoa gmp boost_filesystem boost_system boost_program_options)
link_directories (/usr/local/lib)
include_directories(/usr/local/include/wedge-0.3)

add_executable(demonblast ${SOURCES})
add_executable(debug EXCLUDE_FROM_ALL ${SOURCES} )
target_compile_options(demonblast PUBLIC -O3)
target_compile_options(debug PUBLIC -O0)

#add targets in subdirectory programs, if it exists. This can be used to build other programs that use the code of demonblast
if (EXISTS programs)
	add_subdirectory(programs EXCLUDE_FROM_ALL)
endif()
add_subdirectory(test)



