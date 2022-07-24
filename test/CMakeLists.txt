cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
add_compile_options(-O0)

set (TESTPROGRAMS test_automorphisms test_gauss test_linear_solve test_listofarrows test_matrixbuilder test_niceliegroup test_partitions test_permutations test_signs test_tablemode test_tree)
enable_testing()

foreach(test ${TESTPROGRAMS})
	add_executable(${test} ${test}.cpp ${CMAKE_SOURCE_DIR}/log.cpp)
	add_test( NAME run${test} COMMAND ${test} )
endforeach()

file(GLOB correct_outputs ${CMAKE_SOURCE_DIR}/test/*.ok)


foreach (output ${correct_outputs})
	message(${output})
	get_filename_component(output_base ${output} NAME_WLE)
	add_test( NAME ${output} COMMAND ${CMAKE_COMMAND} -E compare_files  ${CMAKE_CURRENT_BINARY_DIR}/${output_base}.out  ${CMAKE_SOURCE_DIR}/test/${output_base}.ok )

endforeach()
