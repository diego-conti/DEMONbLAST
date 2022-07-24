#include "gauss.cpp"
#include "matrixbuilder.h"
#include <cassert>
#include <iostream>

using namespace std;
using namespace GiNaC;

void test_matrix(Matrix m) {
	assert(m.rows()==2);
	assert(m.cols()==3);
	for (int i=0;i<m.rows();++i)
	for (int j=0;j<m.cols();++j)
		assert(m(i,j)==4);
}

void test_ginac_matrix() {
	auto m=matrix{2,2,lst{4,4,4,4}};
	auto m2=matrix{2,1,lst{4,4}};
	test_matrix(adjoin(MatrixBuilder{m},MatrixBuilder{m2}));
}

void test_transpose_matrix() {
	auto m=matrix{2,2,lst{4,4,4,4}};
	auto m2=matrix{1,2,lst{4,4}};
	test_matrix(adjoin(transpose(MatrixBuilder{m}),transpose(MatrixBuilder{m2})));
}

void test_constant_matrix_builder2() {
	auto m=adjoin(ConstantMatrixBuilder{2,2,4},ConstantMatrixBuilder{2,1,4});
	test_matrix(m);
}

void test_constant_matrix_builder() {
	auto m=adjoin(ConstantMatrixBuilder{2,3,4});
	test_matrix(m);
}

int main() {
  cout<<"testing constant matrix builder...";
  test_constant_matrix_builder();
  test_constant_matrix_builder2();
  cout<<"OK"<<endl;
  cout<<"testing ginac matrix builder...";
  test_ginac_matrix();
  cout<<"OK"<<endl;
  cout<<"testing transpose matrix builder...";
  test_transpose_matrix();
  cout<<"OK"<<endl;
}
