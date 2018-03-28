#include "../gauss.cpp"
#include <cassert>
#include <iostream>
#include <wedge/ginaclinalg.h>

using namespace std;
using namespace GiNaC;
using namespace Wedge;


std::default_random_engine generator;
std::uniform_int_distribution<int> distribution(-100,100);

void fill_at_random(Matrix& m) {
  for (int i=0;i<m.rows();++i)
  for (int j=0;j<m.cols();++j)
    m(i,j)=distribution(generator);
}

void fill_at_random(GinacLinAlgAlgorithms::IndependenceMatrix& m) {
  for (int i=0;i<m.rows();++i)
  for (int j=0;j<m.cols();++j)
    m.M(i,j)=distribution(generator);
}
 
vector<int>  independent_rows_over_Q(GinacLinAlgAlgorithms::IndependenceMatrix& M) {
  M.ChooseLinearlyIndependentRows();
  return {M.IndependentRowsBegin(),M.IndependentRowsEnd()};
}

//using MatrixType = Matrix;
using MatrixType = GinacLinAlgAlgorithms::IndependenceMatrix;

int main() {
  for (int n=3;n<=30;++n)
  for (int m=3;m<=12;++m)
  for (int k=0;k<100;++k) {
    MatrixType M(n,m);
    fill_at_random(M);
    cout<<independent_rows_over_Q(M)<<endl;
   }
}
