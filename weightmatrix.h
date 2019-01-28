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
	bool operator<(Z2 other) const {return x==0 && other.x==1;}
	operator ex() const {return x;}
	int to_Z_star() const {return x? -1 : 1;}
};

inline ostream& operator<<(ostream& os, Z2 x) {
	return os<< (x==Z2{}? 0 : 1);
}

vector<Z2> sign_configuration_to_vector(int dimension, const SignConfiguration& epsilon);


class WeightIterator {
	const vector<WeightAndCoefficient>& weights;
	vector<int>::const_iterator i;
public:
	WeightIterator(const vector<WeightAndCoefficient>& weights,	vector<int>::const_iterator i) : weights{weights},i{i} {}
	const WeightAndCoefficient& operator*() const {return weights[*i];}
	const WeightAndCoefficient* operator->() const {return &weights[*i];}
	bool operator!=(const WeightIterator& j) const {return i!=j.i;}
	bool operator==(const WeightIterator& j) const {return i==j.i;}
	WeightIterator& operator++() {++i; return *this;}
	WeightIterator operator++(int) {
		WeightIterator pre{weights,i};
		++i;
		return pre;
	}
};


//The matrix M_Delta, called "root matrix" in the literature
//Computes internally a maximal set of linearly independent rows to mark weights where parameters or signs can be eliminated, without affecting order
//internally, rows a_1,...,a_n ordered in such a way that the last k rows are basis of Span{a_1,..,a_n} over Z_2 and the last h rows span are a basis  of Span{a_1,..,a_n} over Q
class WeightMatrix {
	friend class SignConfigurations;
	int independent_rows_over_Z2=0, independent_rows_over_Q=0;
	vector<WeightAndCoefficient> weights;
	int cols_;
	vector<int> basis_over_Z2, complement_of_basis_over_Z2;
  int sigma_on_VDelta(const vector<int>& sigma, const WeightAndCoefficient& weight) const {
  	Weight sigmaw{sigma[weight.node_in1], sigma[weight.node_in2],sigma[weight.node_out]};
  	if (sigmaw.node_in1>sigmaw.node_in2) swap(sigmaw.node_in1,sigmaw.node_in2);
  	auto it=find_if(weights.begin(),weights.end(),[sigmaw](auto& w) {return w.same_weight_as(sigmaw);});
  	assert (it!=weights.end());
  	return it-weights.begin();
	}
  Matrix submatrix_JDelta2() const {
  	Matrix submatrix{independent_rows_over_Z2,cols_};
	 	auto weight=Z2basis_begin();
 	  for (int i=0;i<submatrix.rows();++i)
 			populate_row(i, *weight++,submatrix);
 		assert (weight==Z2basis_end());
 		return submatrix;
  }
  Matrix submatrix_IDelta_setminus_JDelta2() const {
    Matrix submatrix{weights.size()-independent_rows_over_Z2,cols_};
	 	auto weight=complement_of_Z2basis_begin();
 	  for (int i=0;i<submatrix.rows();++i)
 			populate_row(i, *weight++,submatrix);
 		assert (weight==complement_of_Z2basis_end());
 		return submatrix;
	}
public:
	WeightMatrix(vector<WeightAndCoefficient>&& weights,int dimension);
	template<typename Container> 
	WeightMatrix(Container&& weights,int dimension) : WeightMatrix{{forward<Container>(weights).begin(),forward<Container>(weights).end()},dimension} {}
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
  vector<int> sigma_on_VDelta(const vector<int>& sigma) const {
		vector<int> result;
		transform(weights.begin(),weights.end(),back_inserter(result),
			[this,&sigma]	(auto& weight) {return this->sigma_on_VDelta(sigma,weight);}
		);
		return result;
  }
  WeightIterator Z2basis_begin() const {return WeightIterator{weights,basis_over_Z2.begin()};}
  WeightIterator Z2basis_end() const {return WeightIterator{weights,basis_over_Z2.end()};}
  WeightIterator complement_of_Z2basis_begin() const {return WeightIterator{weights,complement_of_basis_over_Z2.begin()};}
  WeightIterator complement_of_Z2basis_end() const {return WeightIterator{weights,complement_of_basis_over_Z2.end()};}
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
		transform(MDelta.Z2basis_begin(),MDelta.Z2basis_end(),result.begin(),
			[this,&sigma,&w] (Weight weight) {
				auto sign_and_index=sigma_weight(sigma,weight);
				return sign_and_index.first+w[sign_and_index.second];
			});
		nice_log<<"sigma_w_projection_to_JDelta2= "<<horizontal(result)<<endl;
		return result;
  }
  vector<Z2> sigma_w_projection_to_IDelta_minus_JDelta2(const vector<int>& sigma, const vector<Z2>& w) {
  	vector<Z2>  result(MDelta.rows()-MDelta.rank_over_Z2());
		transform(MDelta.complement_of_Z2basis_begin(),MDelta.complement_of_Z2basis_end(),result.begin(),
			[this,&sigma,&w] (Weight weight) {
				auto sign_and_index=sigma_weight(sigma,weight);
				return sign_and_index.first+w[sign_and_index.second];
			});
		nice_log<<"sigma_w_projection_to_IDelta_minus_JDelta2= "<<horizontal(result)<<endl;
		return result;
  }
  
  SignConfiguration delta(const vector<int>& sigma, const SignConfiguration& epsilon) {
  		auto w_epsilon=sign_configuration_to_vector(MDelta.rows(),epsilon);
  		auto sigma_w_projection_to_JDelta2=this->sigma_w_projection_to_JDelta2(sigma,w_epsilon);
  		auto variables=	 generate_variables<Unknown>(N.a,MDelta.cols());
	 		auto x=solve_over_Z2(
	 			complete_matrix(MDelta.submatrix_JDelta2(),exvector{sigma_w_projection_to_JDelta2.begin(),sigma_w_projection_to_JDelta2.end()}),
				variables
	 		);
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
			return result;
  	} 
public:
	SignConfigurations(const WeightMatrix& MDelta,const list<vector<int>>& nontrivial_automorphisms) :
		 MDelta{MDelta}, sign_configurations_{SignConfiguration::all_configurations(MDelta.rows()-MDelta.rank_over_Z2())} {
		factor_out_automorphisms(nontrivial_automorphisms);
	}
	
	list<SignConfiguration> sign_configurations() && {return move(sign_configurations_);}
};

/*
class ImageMod2 {
	map<vector<Z2>,list<vector<Z2>>> preimages;
	static string to_string(const vector<Z2>& vec, const list<vector<Z2>>& preimages) {
		string result=horizontal(vec)+" from:\n";
		for (auto& v : preimages) result+="\t"+horizontal(v)+"\n";
		return result;
	}
public:
	string to_string() const {
		string result;
		for (auto& vector_and_preimages : preimages) result+=to_string(vector_and_preimages.first, vector_and_preimages.second)+"\n";
		return result;
	}
	void insert(const vector<Z2>& delta, const vector<Z2>& imdelta) {
		preimages[imdelta].push_back(delta);
	}
};
*/

class ImageMod2 {
	map<vector<Z2>,list<string>> preimages;
	static string to_string(const vector<Z2>& vec, const list<string>& preimages) {
		string result=horizontal(vec)+" from:\n";
		for (auto& v : preimages) result+="\t"+v+"\n";
		return result;
	}
public:
	string to_string() const {
		string result;
		for (auto& vector_and_preimages : preimages) result+=to_string(vector_and_preimages.first, vector_and_preimages.second);
		return result;
	}
	void insert(const string& delta, const vector<Z2>& imdelta) {
		preimages[imdelta].push_back(delta);
	}
};

ImageMod2 image_mod2(const WeightMatrix& weight_matrix);

#endif
