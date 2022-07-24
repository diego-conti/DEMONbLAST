#include "gauss.cpp"
#include <cassert>
#include <iostream>
#include <ginac/ginac.h>

using namespace std;
using namespace GiNaC;

void test_gauss() {
  assert ( (independent_rows_over_Q(Matrix{{1,2,3},{2,3,1},{3,1,2}}) == vector<int>{0,1,2} ));
  assert ( (independent_rows_over_Z2(Matrix{{1,2,3},{2,3,1},{3,1,2}}) == vector<int>{0,1} ));
  assert ( (independent_rows_over_Q(Matrix{{1,2,3},{1,2,3},{2,3,1},{3,1,2}}) == vector<int>{0,2,3} ));
  assert ( (independent_rows_over_Z2(Matrix{{1,2,3},{1,4,3},{2,3,1},{3,1,2}}) == vector<int>{0,2} ));
}

void test_matrix() {
  Matrix m{{1,0,0},{0,1,0},{0,0,1}};
  assert(m(0,0)==1);
  assert(m(1,0)==0);
  assert(m(2,0)==0);
  assert(m(2,2)==1);
  assert(m.rows()==3);
  assert(m.cols()==3);
  m(0,0)=2;
  assert(m(0,0)==2);
 }


exvector three_variables() {
  return {symbol("x"),symbol("y"),symbol("z")};
}
exvector variables(int n) {
  exvector result;
  for (int i=0;i<n;++i) result.push_back(realsymbol("n"+to_string(i)));
  return result;
}


matrix to_ginac_matrix(const Matrix& matrix) {
  GiNaC::matrix m(matrix.rows(),matrix.cols());
  for (int i=0;i<matrix.rows();++i)
  for (int j=0;j<matrix.cols();++j)
    m(i,j)=matrix(i,j);
   return m;
}

matrix to_ginac_matrix(const exvector& v) {
  matrix m(v.size(),1);
  for (int i=0;i<m.rows();++i)
    m(i,0)=v[i];
   return m;
}

bool is_even(ex x) {  
  x=x.expand().subs(wild()*2==0, subs_options::algebraic);
  if (!is_a<numeric>(x)) return false;
  return ex_to<numeric>(x).to_int()%2==0;  
}

void test_solve_over_Z2(const Matrix& complete_matrix) {
  int no_variables=complete_matrix.cols()-1;
  auto sol=solve_over_Z2(complete_matrix,variables(no_variables));
  sol.push_back(0);
  auto ginac_matrix= to_ginac_matrix(complete_matrix);
  auto ginac_vector=to_ginac_matrix(sol);
  auto product=ex_to<matrix>((ginac_matrix*ginac_vector).evalm());
  for (int i=0;i<ginac_matrix.rows();++i) {
    assert(is_even(ginac_matrix(i,no_variables)-product(i,0)));
}
}


void test_solve_over_Z2() {
  assert(solve_over_Z2({{1,2,3,0},{4,5,6,0},{7,8,9,1}},three_variables()).empty());  
  assert(!solve_over_Z2({{1,2,3,0},{4,5,6,0},{7,8,9,0}},three_variables()).empty());  
  auto sol=solve_over_Z2({{1,2,3,1},{4,5,6,0},{7,8,9,1}},three_variables());
  assert(ex_to<numeric>(sol[0]+sol[2]).to_int()%2==1);
  assert(ex_to<numeric>(sol[1]).to_int()%2==0);  
  test_solve_over_Z2(
    {{1, 0, -1, 0, 0, 0, 0, -1, 0, 1},
     {1, 0, 0, -1, 0, 0, 0, 0, -1,   0},
      {1, 0, 0, 0, -1, 0, -1, 0, 0, 1},
      {0, 1, 0, 0, -1, 0, 0, 0, -1,   1},
      {0, 0, 1, 0, 0, -1, 0, -1, 0, 1},
       {0,  0, 0, 1, 0, -1, 0,   0, -1, 0},
        {0, 1, 0, 0, 0, 0, -1, -1, 0, 0},
         {0, 0, 0, 0, 1,   0, -1, 0, -1, 1},
          {0, 0, 0, 0, 0, 1, 0, -1, -1, 0} }
   );
  test_solve_over_Z2(
    {{1,0,-1,0,0,0,0,-1,0,1},
    {1,0,0,-1,0,0,0,0,-1,0},
    {1,0,0,0,-1,0,-1,0,0,1},
    {0,1,0,0,-1,0,0,0,-1,1},
    {0,0,1,0,0,-1,0,-1,0,1},
    {0,0,0,1,0,-1,0,0,-1,0},
    {0,1,0,0,0,0,-1,-1,0,0},
    {0,0,0,0,1,0,-1,0,-1,1},
    {0,0,0,0,0,1,0,-1,-1,0}
  });


}


int main() {
  cout<<"testing matrix...";
  test_matrix();
  cout<<"OK"<<endl;
  
  cout<<"testing gauss...";
  test_gauss();
  cout<<"OK"<<endl;
  
  cout<<"testing solve_over_Z2...";
  test_solve_over_Z2();
  cout<<"OK"<<endl;
}
