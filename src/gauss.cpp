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
#include "gauss.h"
#include "log.h"

using namespace std;
using namespace GiNaC;

thread_local ex ex_0{0};
Matrix::Matrix(int rows, int cols) :  no_rows{rows}, no_cols{cols}, flat_representation(rows*cols, ex_0) {}


Matrix::Matrix(std::initializer_list<std::initializer_list<ex>> elements) {
    no_rows=elements.size();    
    assert(no_rows>0);
    for (auto row : elements)
      flat_representation.insert(flat_representation.end(),row.begin(),row.end());
    no_cols=flat_representation.size()/no_rows;
}

Matrix::Matrix(const matrix& m) : Matrix(m.rows(),m.cols()) {
  for (int i=0;i<rows();++i)
  for (int j=0;j<cols();++j)
    operator()(i,j)=m(i,j);
}


void print_row(ostream& os, const Matrix& m, int row) {
    os<<"{"<<m(row,0);
    int j=0;
    while (++j<m.cols()) os<<','<<m(row,j);
    os<<"}"<<endl;  
}

ostream& operator<<(ostream& os, const Matrix& m) {
  os<<"{";
  for (int i=0;i<m.rows();++i)
    print_row(os,m,i);
  return os<<"}"<<endl;  
}

exvector Matrix::column(int j) const {
  exvector result;
  for (int i=0;i<rows();++i) result.push_back((*this)(i,j));
  return result;
}

exvector Matrix::image_of(const exvector& coefficients) const {
  exvector result(rows());
  for (int i=0;i<rows();++i)
  for (int j=0;j<cols();++j)
    result[i]+=coefficients[j]*(*this)(i,j);
  return result;
}

class PivotableMatrix : private Matrix {
  friend class LinearEquation;
  friend std::ostream& operator<<(std::ostream& os, const PivotableMatrix& m);
public:  
  using PivotPosition = pair<int,int>;

  template<class Reducer, class ... Other>
  PivotableMatrix(Matrix&& matrix, Reducer&& reducer, Other... other) : Matrix{move(matrix)}  {
    init_actual_column_position();
    reduce(forward<Reducer>(reducer), other...);
  }
  template<class Reducer, class ... Other>
  PivotableMatrix(const Matrix& matrix, Reducer&& reducer, Other... other) : Matrix{matrix} {
    init_actual_column_position();
    reduce(forward<Reducer>(reducer),other...);
  }
  PivotableMatrix(int rows, int cols) : Matrix{rows,cols} {
      init_actual_column_position();
  }

  template<typename Reducer>      
  void eliminate_past_row(int pivot_row, Reducer&& reducer) {
      int pivot_col=pivots_used++;
      ex pivot=at(pivot_row,pivot_col);      
			for (unsigned i=pivot_row+1; i<rows(); ++i) {
				for (unsigned j=pivot_col+1; j<cols(); ++j)
					at(i,j)=reducer((pivot*at(i,j) - at(i,pivot_col)*at(pivot_row,j))/last_pivot);
			  at(i,pivot_col)=0;
			}
      last_pivot=pivot;
  }
  
  bool pivot_in_row(int i) {
    int j=first_non_zero_element_in_row(i);
    if (j==cols()) return false;
  	swap_columns(j,pivots_used);
    return true;  	
  }
  bool pivot_in_row(int i, ignore_last_column_t) {
    int j=first_non_zero_element_in_row(i);
    if (j>=cols()-1) return false;
  	swap_columns(j,pivots_used);
    return true;  	
  }  
  
  using Matrix::rows;
  using Matrix::cols;
  
  template<typename... Flags>  vector<PivotPosition> pivots(Flags... flags) const {
    vector<PivotPosition> result;
    int row=0;
    for (int pivot=0;pivot<pivots_used;++pivot) {
      while (row_is_zero(row, flags...)) ++row;
      result.emplace_back(row,actual_column_position[pivot]);
      ++row;
    }
    return result;
  }
  
 exvector last_column() const  {
  exvector last_column;
  int last_col=cols()-1;
  for (int i=0;i<rows();++i)
    last_column.push_back(at(i,last_col));
  return last_column;
  }  

private:
  bool row_is_zero(int row) const {
    for (int i=0;i<cols();++i) if (!Matrix::operator()(row,i).is_zero()) return false;
    return true;
  }
  bool row_is_zero(int row,ignore_last_column_t) const {
    for (int i=0;i<cols()-1;++i) if (!Matrix::operator()(row,i).is_zero()) return false;
    return true;
  }
  int first_non_zero_element_in_row(int i) const {
  	  int k = pivots_used;
	  	while (k<cols() && at(i,k).is_zero())
  			++k;
  		return k;
  }
  void swap_columns(int c1,int c2) {
    if (c1!=c2) swap(actual_column_position[c1],actual_column_position[c2]);
  }
  ex& at(int zero_based_index_i, int zero_based_index_j) {return Matrix::operator() (zero_based_index_i,actual_column_position[zero_based_index_j]);}
  ex at(int zero_based_index_i, int zero_based_index_j) const {return Matrix::operator() (zero_based_index_i,actual_column_position[zero_based_index_j]);}
  void init_actual_column_position() {
      actual_column_position.resize(cols());
			for (int i=0;i<cols();++i) actual_column_position[i]=i;  
  }
  ex last_pivot=1;
  int pivots_used=0;
  vector<int> actual_column_position;
};

std::ostream& operator<<(std::ostream& os, const PivotableMatrix& m) {
  os<<"without reordering columns "<<endl<<static_cast<const Matrix&>(m);
  os<<"with reordering"<<endl;
  for (int i=0;i<m.rows();++i) {
  for (int j=0;j<m.cols();++j)
    os<<m.at(i,j)<<",";
  os<<endl;
  }
  return os;
}


template<typename Reducer>
vector<int> independent_rows_in(Matrix&& m, Reducer&& reducer) {
  vector<int> independent_rows;
  independent_rows.reserve(m.rows());
  PivotableMatrix M{move(m),reducer};
	for (unsigned i=0; i<M.rows(); ++i)
		if (M.pivot_in_row(i)) {
			independent_rows.push_back(i);
			M.eliminate_past_row(i,reducer);	
   }
  return independent_rows;
}

vector<int> independent_rows_over_Q(Matrix&& m) {
  static auto reducer = [] (ex x) {return x.expand();};
  return independent_rows_in(move(m),reducer);
}
vector<int> independent_rows_over_Q(const Matrix& m) {
  return independent_rows_over_Q(Matrix{m});
}

ex reduce_mod_Z2(ex x) {
  assert(is_a<numeric>(x));
  const numeric& as_numeric = ex_to<numeric>(x);
//  assert(as_numeric.is_integer());
  return as_numeric.to_int()%2;
}


vector<int> independent_rows_over_Z2(Matrix&& m) {
  return independent_rows_in(move(m),reduce_mod_Z2);
}

vector<int> independent_rows_over_Z2(const Matrix& m) {
  return independent_rows_over_Z2(Matrix{m});
}


template<typename Reducer>
class EchelonReducedMatrix : public Matrix {

public:

};

class LinearEquation {
  const PivotableMatrix& matrix;
  int row;
  int no_variables;
public:
  LinearEquation(const PivotableMatrix& complete_matrix,int row) : matrix{complete_matrix}, row{row}, no_variables{matrix.cols()-1} {}
  ex solve_for(int column_with_pivot, const exvector& partial_solution) const {
    ex pivot=matrix(row,column_with_pivot);
    if (pivot.is_zero()) {
      nice_log<<matrix<<endl;
      nice_log<<row<<endl;
   }
    assert (!pivot.is_zero());
    assert(partial_solution[column_with_pivot].is_zero());
    ex sum=evaluate_on(partial_solution);
    return -sum/pivot;
  }
private:
  ex evaluate_on(const exvector& variables) const {
    assert(no_variables==variables.size());
    ex result;
    for (int col=0;col<no_variables;++col)
      result+=variables[col]*matrix(row,col);
    return result-matrix(row,no_variables);
  }
};


bool has_solution(const PivotableMatrix& reduced_complete_matrix) {
  exvector last_column=reduced_complete_matrix.last_column();
 for (auto pivot :  reduced_complete_matrix.pivots(Matrix::ignore_last_column)) {
    int row=pivot.first;
    last_column[row]=0;
  }
  return all_of(last_column.begin(),last_column.end(),[](ex x) {return x.is_zero();});
}

exvector partial_solution(const PivotableMatrix& reduced_complete_matrix, const exvector& variables) {
  assert(variables.size()==reduced_complete_matrix.cols()-1);
  exvector partial_solution = variables;
  for (auto pivot :  reduced_complete_matrix.pivots(Matrix::ignore_last_column)) {
    int column=pivot.second;
    partial_solution[column]=0;
  }
  return partial_solution;
}

exvector solve_by_back_substitution(const PivotableMatrix& reduced_complete_matrix, const exvector& variables) {
  exvector partial_solution = ::partial_solution(reduced_complete_matrix,variables);
  auto pivots = reduced_complete_matrix.pivots(Matrix::ignore_last_column);
  for (auto i =pivots.rbegin();i!=pivots.rend();++i ) {  
    int pivot_row=i->first, pivot_column=i->second;
    LinearEquation eq(reduced_complete_matrix,pivot_row);
    partial_solution[pivot_column]=eq.solve_for(pivot_column,partial_solution);
  }
  return partial_solution;
} 



template<typename Reducer> 
exvector solve(const Matrix& complete_matrix, const exvector& variables, Reducer&& reducer) {
  PivotableMatrix echelon_reduced_complete_matrix{complete_matrix,reducer,Matrix::ignore_last_column};
	for (unsigned i=0; i<echelon_reduced_complete_matrix.rows(); ++i)
		if (echelon_reduced_complete_matrix.pivot_in_row(i,Matrix::ignore_last_column))
			echelon_reduced_complete_matrix.eliminate_past_row(i,reducer);
  for (auto x : echelon_reduced_complete_matrix.pivots(Matrix::ignore_last_column)) 
	    nice_log<<"pivot "<<x.first<<" "<<x.second<<endl;	  
	if (!has_solution(echelon_reduced_complete_matrix)) {
	  nice_log<<"no solution to"<<endl<<echelon_reduced_complete_matrix<<endl;
	  for (auto x : echelon_reduced_complete_matrix.pivots(Matrix::ignore_last_column)) 
	    nice_log<<"pivot "<<x.first<<" "<<x.second<<endl;	  
	  return {};
	 }
  return solve_by_back_substitution(echelon_reduced_complete_matrix,variables);
}

exvector solve_over_Q(const Matrix& complete_matrix, const exvector& variables) {
  static auto reducer = [] (ex x) {return x.expand();};
  return solve(complete_matrix,variables, reducer);
}

exvector solve_over_Z2(const Matrix& complete_matrix, const exvector& variables) {
  auto solution = solve(complete_matrix,variables, reduce_mod_Z2);
  for (auto& variable : solution) {
    variable=variable.expand().subs(wild()*2==0, subs_options::algebraic);
    if (is_a<numeric>(variable)) variable=reduce_mod_Z2(variable);
  }
  return solution;
}

Matrix complete_matrix(const Matrix& incomplete_matrix, exvector constant_terms) {
   assert(incomplete_matrix.rows()==constant_terms.size());
    Matrix M(incomplete_matrix.rows(),incomplete_matrix.cols()+1);
		for (int i=0;i<incomplete_matrix.rows();++i)
    for (int j=0;j<incomplete_matrix.cols();++j) {
       M(i,j)=incomplete_matrix(i,j);
       M(i,incomplete_matrix.cols())=constant_terms[i];
    }
	  return M;	  
}

