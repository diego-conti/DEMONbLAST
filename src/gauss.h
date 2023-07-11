/*  Copyright (C) 2018-2023 by Diego Conti, diego.conti@unipi.it      
                                                                     
    This file is part of DEMONbLAST
	                                                                     
    DEMONbLAST is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    DEMONbLAST is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with DEMONbLAST.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef GAUSS_H
#define GAUSS_H


#include "includes.h"


using std::vector;
using GiNaC::ex;
using GiNaC::exvector;
using GiNaC::matrix;

class Matrix {
public:
  Matrix(int rows, int cols);
  Matrix(std::initializer_list<std::initializer_list<ex>> elements);
  explicit Matrix(const matrix& m);
  ex at (int zero_based_index_i, int zero_based_index_j) const {return flat_representation[zero_based_index_i*no_cols+zero_based_index_j];}
  ex operator() (int zero_based_index_i, int zero_based_index_j) const {return flat_representation[zero_based_index_i*no_cols+zero_based_index_j];}
  ex& operator() (int zero_based_index_i, int zero_based_index_j) {return flat_representation[zero_based_index_i*no_cols+zero_based_index_j];}
  int rows() const {return no_rows;}
  int cols() const {return no_cols;}
  class ignore_last_column_t {};
  constexpr static ignore_last_column_t ignore_last_column{};
  
  exvector column(int j) const;
  exvector row(int i) const;
  exvector image_of(const exvector& coefficients) const;
protected:
  template<typename Reducer>
  void reduce(Reducer&& reducer) {for (auto& x : flat_representation) x=reducer(x);}
  template<typename Reducer>
  void reduce(Reducer&& reducer, ignore_last_column_t) {
    for (int i=0;i<no_rows;++i) 
    for (int j=0;j<no_cols-1;++j) 
      operator() (i,j)=reducer(operator()(i,j));
  }
private:
  int no_rows, no_cols;
  exvector flat_representation;
};

vector<int> independent_rows_over_Q(const Matrix& m);
vector<int> independent_rows_over_Q(Matrix&& m);
vector<int> independent_rows_over_Z2(const Matrix& m);
vector<int> independent_rows_over_Z2(Matrix&& m);
exvector solve_over_Z2(const Matrix& complete_matrix, const exvector& variables);
exvector solve_over_Q(const Matrix& complete_matrix, const exvector& variables);

std::ostream& operator<<(std::ostream& os, const Matrix& m);

ex reduce_mod_Z2(ex x);

Matrix complete_matrix(const Matrix& incomplete_matrix, exvector constant_terms);
template<typename Parameter, typename NameClass> exvector generate_variables(const NameClass& n, int max_index) {
   exvector result;
   for (int i=1;i<=max_index;++i) result.emplace_back(Parameter{n(i)});
  return result;
}

#endif
