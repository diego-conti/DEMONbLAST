/*  Copyright (C) 2018 by Diego Conti, diego.conti@unimib.it      
                                                                     
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
#ifndef WEIGHTMATRIX_H
#define WEIGHTMATRIX_H

#include "matrixbuilder.h"
#include "weightbasis.h"

template<typename Matrix>
void populate_row(int row, Weight weight, Matrix& matrix) {
     ex ex_1=ex{1};
     ex ex_minus_1=ex{-1};
		matrix(row,weight.node_out)=ex_1;
 		matrix(row,weight.node_in1)=ex_minus_1;
   	matrix(row,weight.node_in2)=ex_minus_1;
}
	
WEDGE_DECLARE_NAMED_ALGEBRAIC(Unknown,realsymbol);

class Z2 {
	int x;
public:
	Z2(int x=0) : x{x%2} {}
	bool operator==(Z2 other) const {return x==other.x;}
	bool operator!=(Z2 other) const {return x==other.x;}
	Z2 operator+ (Z2 other) const {return {x+other.x};}	
	Z2& operator=(int other) {x=other%2; return *this;}
	Z2& operator=(const Z2& other) =default;
	operator ex() const {return x;}
	int to_Z_star() const {return x? -1 : 1;}
};

inline ostream& operator<<(ostream& os, Z2 x) {
	return os<< (x==Z2{}? 0 : 1);
}

vector<Z2> sign_configuration_to_vector(int dimension, const SignConfiguration& epsilon);

//The matrix M_Delta, called "root matrix" in the literature, with rows a_1,...,a_n ordered in such a way that the last k rows are basis of Span{a_1,..,a_n} over Z_2 and the last h rows span are a basis  of Span{a_1,..,a_n} over Q
class WeightMatrix {
	int independent_rows_over_Z2=0, independent_rows_over_Q=0;
	vector<WeightAndCoefficient> weights;
	int cols_;
public:
	WeightMatrix(vector<WeightAndCoefficient> unordered_weights,int dimension);
	template<typename Container> 
	WeightMatrix(Container&& unordered_weights,int dimension) : WeightMatrix{{unordered_weights.begin(),unordered_weights.end()},dimension} {}
  int rank_over_Z2() const {
		return independent_rows_over_Z2;
  }
  int rank_over_Q() const {
		return independent_rows_over_Q;
  }
  Matrix M_Delta() const {
  	Matrix matrix{rows(),cols()};
	 	auto weight=weights.begin();
 	  for (int i=0;i<matrix.rows();++i)
 			populate_row(i, *weight++,matrix);
 		return matrix;  	
 	}
  Matrix submatrix_JDelta2() const {
  	Matrix submatrix{independent_rows_over_Z2,cols_};
	 	auto weight=weights.begin()+weights.size()-independent_rows_over_Z2;
 	  for (int i=0;i<submatrix.rows();++i)
 			populate_row(i, *weight++,submatrix);
 		return submatrix;
  }
  Matrix submatrix_IDelta_setminus_JDelta2() const {
    Matrix submatrix{weights.size()-independent_rows_over_Z2,cols_};
	 	auto weight=weights.begin();
 	  for (int i=0;i<submatrix.rows();++i)
 			populate_row(i, *weight++,submatrix);
 		return submatrix;
	}
  auto weight_begin() const {return weights.begin();}
  auto weight_end() const {return weights.end();} 	
  int rows() const {return weights.size();}
  int cols() const {return cols_;}
  vector<Z2> image_of(const vector<Z2> v) const {
  	vector<Z2> result;
  	result.reserve(weights.size());
  	for (auto x: weights) result.push_back(v[x.node_in1] + v[x.node_in2]+v[x.node_out]);
  	return result;
  }
};

exvector X_solving_nonRicciflat_Einstein(const WeightMatrix& MDelta);
exvector X_solving_Ricciflat(const WeightMatrix& MDelta);
exvector X_solving_nilsoliton(const WeightMatrix& MDelta);

class SignConfigurations {
 	const WeightMatrix& MDelta;
 	list<SignConfiguration> sign_configurations_;
 	
  void factor_out_automorphisms(const list<vector<int>>& nontrivial_automorphisms) {
  	if (sign_configurations_.empty()) return;
  	for (auto epsilon = sign_configurations_.begin();epsilon!=sign_configurations_.end();++epsilon) {
  		for (auto& sigma : nontrivial_automorphisms) {
  			auto delta_sigma_epsilon=delta(sigma, *epsilon);
  			if (delta_sigma_epsilon!=*epsilon) 
	  			sign_configurations_.remove(delta_sigma_epsilon);
	  	}
  	}
  }
  
  pair<Z2,int> sigma_weight(const vector<int>& sigma,Weight weight) const {
  	Z2 logsign=0;
  	weight={sigma[weight.node_in1],sigma[weight.node_in2],sigma[weight.node_out]};
  	if (weight.node_in1>weight.node_in2) {
  		logsign=1; 
  		swap(weight.node_in1,weight.node_in2);
  	}
  	int index= find_if(MDelta.weight_begin(),MDelta.weight_end(),
  			[weight] (auto& w) {return w.same_weight_as(weight);}
			)-MDelta.weight_begin();
  	assert(index>=0 && index<=MDelta.rows());
  	return {logsign,index};
  }
  	
	vector<Z2> sigma_w_projection_to_JDelta2(const vector<int>& sigma, const vector<Z2>& w) {
		vector<Z2> result(MDelta.rank_over_Z2());
		auto Z2basis_begin=MDelta.weight_begin()+ MDelta.rows()-MDelta.rank_over_Z2();
		transform(Z2basis_begin,MDelta.weight_end(),result.begin(),
			[this,&sigma,&w] (Weight weight) {
				auto sign_and_index=sigma_weight(sigma,weight);
				return sign_and_index.first+w[sign_and_index.second];
			});
		nice_log<<"sigma_w_projection_to_JDelta2= "<<horizontal(result)<<endl;
		return result;
  }
  vector<Z2> sigma_w_projection_to_IDelta_minus_JDelta2(const vector<int>& sigma, const vector<Z2>& w) {
  	vector<Z2>  result(MDelta.rows()-MDelta.rank_over_Z2());
		auto Z2basis_begin=MDelta.weight_begin()+ MDelta.rows()-MDelta.rank_over_Z2();
		transform(MDelta.weight_begin(),Z2basis_begin,result.begin(),
			[this,&sigma,&w] (Weight weight) {
				auto sign_and_index=sigma_weight(sigma,weight);
				return sign_and_index.first+w[sign_and_index.second];
			});
		nice_log<<"sigma_w_projection_to_IDelta_minus_JDelta2= "<<horizontal(result)<<endl;
		return result;
  }
  
  SignConfiguration delta(const vector<int>& sigma, const SignConfiguration& epsilon) {
  		auto w_epsilon=sign_configuration_to_vector(MDelta.rows(),epsilon);
  		nice_log<<"sigma = "<<horizontal(sigma)<<endl;
  		nice_log<<"epsilon = "<<epsilon<<endl;
  		nice_log<<"w_epsilon = "<<horizontal(w_epsilon)<<endl;
  		auto sigma_w_projection_to_JDelta2=this->sigma_w_projection_to_JDelta2(sigma,w_epsilon);
  		auto variables=	 generate_variables<Unknown>(N.a,MDelta.cols());
	 		auto x=solve_over_Z2(
	 			complete_matrix(MDelta.submatrix_JDelta2(),exvector{sigma_w_projection_to_JDelta2.begin(),sigma_w_projection_to_JDelta2.end()}),
				variables
	 		);
	 		nice_log<<"x = "<<horizontal(x)<<endl;
	 		nice_log<<"M_Delta2,x = "<<horizontal(MDelta.submatrix_JDelta2().image_of(x))<<endl;
	 		auto delta=MDelta.submatrix_IDelta_setminus_JDelta2().image_of(x);	
	 		auto sigma_w=sigma_w_projection_to_IDelta_minus_JDelta2(sigma,w_epsilon);
			SignConfiguration result{delta.size()};
			for (int i=0;i<delta.size();++i) {
				for (auto variable : variables) delta[i]=delta[i].subs(variable==0);
				nice_log<<delta[i]<<endl;
				assert(is_a<numeric>(delta[i]));
				Z2 delta_i{ex_to<numeric>(delta[i]).is_even()? 0 : 1};
				result[i]=(sigma_w[i]+delta_i).to_Z_star();
			}
			nice_log<<"delta(sigma,epsilon)+w_epsilon = "<<result<<endl;
			return result;
  	} 
public:
	SignConfigurations(const WeightMatrix& MDelta,const list<vector<int>>& nontrivial_automorphisms) :
		 MDelta{MDelta}, sign_configurations_{SignConfiguration::all_configurations(MDelta.rows()-MDelta.rank_over_Z2())} {
		factor_out_automorphisms(nontrivial_automorphisms);
	}
	
	list<SignConfiguration> sign_configurations() && {return move(sign_configurations_);}
};


#endif
