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

//The matrix M_Delta, called "root matrix" in the literature
class WeightMatrix {
  matrix weight_matrix;
  vector<int> row_to_parameter_correspondence;
  vector<int> independent_rows_over_Z2;
  vector<int> independent_rows_over_Q;
  vector<WeightAndCoefficient> weights;
  exvector X_ijk_;
  bool derivations_traceless=false;
 	list<SignConfiguration> sign_configurations_;
 	
  Matrix complete_matrix_transposed_M_Delta(ex lambda) const {
    Matrix M(weight_matrix.cols(),weight_matrix.rows()+1);
		for (int i=0;i<weight_matrix.rows();++i)
    for (int j=0;j<weight_matrix.cols();++j) {
       M(j,i)=weight_matrix(i,j);
       M(j,weight_matrix.rows())=lambda;
    }
	  return M;	  
  }
  Matrix reordered_matrix_for_elimination(const vector<WeightAndCoefficient>& weights) {
 	  Matrix M(weight_matrix.rows(),weight_matrix.cols());
 	  int next_row_in_front=0, next_row_in_back=weight_matrix.rows();
		for (int i=0;i<weights.size();++i)
  		if (weights[i].parameter_eliminated()) {
  		  row_to_parameter_correspondence[next_row_in_front]=i;
  		  populate_row(next_row_in_front++,weights[i],M);
  	}
     else {
       	row_to_parameter_correspondence[--next_row_in_back]=i;
  		  populate_row(next_row_in_back,weights[i],M);
  	}
	  return M;
  }

  void store_X_ijk() {
    X_ijk_=solve_over_Q(complete_matrix_transposed_M_Delta(1),generate_variables<StructureConstant>(N.x,weights.size()));
    if (!X_ijk_.empty()) {
      assert(X_ijk_.size()==weights.size());
      for (int i=0;i<weights.size();++i)
        weights[i].X_ijk=X_ijk_[i];
      derivations_traceless=true;
    }  
  }
  int rows_in_JDelta, rows_in_JDelta2;
  void determine_sign_ambiguities() {  
 	 rows_in_JDelta= rows_in_JDelta2=0;
	 	vector<WeightAndCoefficient> reordered_weights;    	
		for (auto row : independent_rows_over_Z2)  {
			weights[row].eliminate_parameter();
			reordered_weights.push_back(weights[row]);
	 		++rows_in_JDelta2;
		}
    independent_rows_over_Q=::independent_rows_over_Q(reordered_matrix_for_elimination(weights));
		for (auto row : independent_rows_over_Q) {
	 	 	++rows_in_JDelta;
		  int parameter=row_to_parameter_correspondence[row];
		  if (!weights[parameter].parameter_eliminated()) {
		  		weights[parameter].eliminate_parameter();
					reordered_weights.push_back(weights[parameter]);					
			}		
		}
		for (auto& weight: weights) 
			if (!weight.parameter_eliminated())
				reordered_weights.push_back(weight);
		weights={reordered_weights.rbegin(),reordered_weights.rend()};	
		sign_configurations_= SignConfiguration::all_configurations(weights.size()-rows_in_JDelta2);
		nice_log<<"rank J_Delta2 = "<<rows_in_JDelta2<<", rank J_Delta="<<rows_in_JDelta<<endl;
  }
  
  pair<Z2,int> sigma_weight(const vector<int>& sigma,Weight weight) const {
  	Z2 logsign=0;
  	weight={sigma[weight.node_in1],sigma[weight.node_in2],sigma[weight.node_out]};
  	if (weight.node_in1>weight.node_in2) {
  		logsign=1; 
  		swap(weight.node_in1,weight.node_in2);
  	}
  	int index= find_if(weights.begin(),weights.end(),
  		[weight] (const WeightAndCoefficient& weight_and_coeff) {return weight_and_coeff.same_weight_as(weight);}
  	)-weights.begin();
  	return {logsign,index};
  }
  
	vector<Z2> sigma_w_projection_to_JDelta2(const vector<int>& sigma, const vector<Z2>& w) {
		vector<Z2> result(rows_in_JDelta2);
		for (int i=0;i<result.size();++i) {
			auto sign_and_index=sigma_weight(sigma,weights[i+weights.size()-rows_in_JDelta2]);
			result[i]=sign_and_index.first+w[sign_and_index.second];
		}
		nice_log<<"sigma_w_projection_to_JDelta2= "<<horizontal(result)<<endl;
		return result;
  }
  vector<Z2> sigma_w_projection_to_IDelta_minus_JDelta2(const vector<int>& sigma, const vector<Z2>& w) {
  	vector<Z2>  result(weights.size()-rows_in_JDelta2);
		for (int i=0;i<result.size();++i) {
			auto sign_and_index=sigma_weight(sigma,weights[i]);
			result[i]=sign_and_index.first+w[sign_and_index.second];
		}
		nice_log<<"sigma_w_projection_to_JDelta_minus_JDelta2= "<<horizontal(result)<<endl;
		return result;
  }
  
  //the submatrix in M_Delta,2 with rows indexed by J_Delta,2
  	Matrix submatrix_JDelta2() const {
  		Matrix submatrix(rows_in_JDelta2,weight_matrix.cols()); 
  		 auto weight=weights.begin()+weights.size()-rows_in_JDelta2;
	 	  for (int i=0;i<submatrix.rows();++i)
	 			populate_row(i, *weight++,submatrix);  	
			nice_log<<"submatrix_JDelta2= "<<submatrix<<endl;
	 		return submatrix;
  	}
  //the submatrix in M_Delta,2 with rows indexed by J_Delta\J_Delta2
  	Matrix submatrix_IDelta_setminus_JDelta2() const {
  		Matrix submatrix(weights.size()-rows_in_JDelta2,weight_matrix.cols()); 
  		auto weight=weights.begin();
	 	  for (int i=0;i<submatrix.rows();++i)
	 			populate_row(i, *weight++,submatrix);  	
			nice_log<<"submatrix_IDelta_setminus_JDelta2= "<<submatrix<<endl;
	 		return submatrix;
  	}
  	 
  	vector<Z2> representative(const SignConfiguration& epsilon) const {
  		vector<Z2> w_epsilon(weights.size());
  		for (int i=0; i<epsilon.size();++i) 
				w_epsilon[i]=epsilon[i]>0? 0 : 1;
  		return w_epsilon;
  	}
  	SignConfiguration delta(const vector<int>& sigma, const SignConfiguration& epsilon) {
  		auto w_epsilon=representative(epsilon);
  		nice_log<<"sigma = "<<horizontal(sigma)<<endl;
  		nice_log<<"epsilon = "<<epsilon<<endl;
  		nice_log<<"w_epsilon = "<<horizontal(w_epsilon)<<endl;
  		auto sigma_w_projection_to_JDelta2=this->sigma_w_projection_to_JDelta2(sigma,w_epsilon);
  		auto variables=	 generate_variables<Unknown>(N.a,weight_matrix.cols());
	 		auto x=solve_over_Z2(
	 			complete_matrix(submatrix_JDelta2(),exvector{sigma_w_projection_to_JDelta2.begin(),sigma_w_projection_to_JDelta2.end()}),
				variables
	 		);
	 		nice_log<<"x = "<<horizontal(x)<<endl;
	 		nice_log<<"M_Delta2,x = "<<horizontal(submatrix_JDelta2().image_of(x))<<endl;
	 		auto delta=submatrix_IDelta_setminus_JDelta2().image_of(x);	
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
  	
  	static matrix matrix_of_ones(int rows) {
  		exvector ones(rows,1);
  		return {rows,1,lst{ones.begin(),ones.end()}};  		
  	}
  	static matrix matrix_of_unknowns(int rows) {
  		auto ones=generate_variables<Unknown>(N.b,rows);
  		return {rows,1,lst{ones.begin(),ones.end()}};  		
  	}
public:
  WeightMatrix(const list<Weight>& weights,int dimension) : 
    weight_matrix{weights.size(),dimension},  
    row_to_parameter_correspondence(weights.size()),
    weights{weights.begin(),weights.end()}
  {
    Matrix matrix_for_elimination{weights.size(),dimension};
    int i=0;
    for (auto weight: weights) {
  		populate_row(i, weight,matrix_for_elimination);    
  		populate_row(i++, weight,weight_matrix);    
  	}
 		independent_rows_over_Z2=::independent_rows_over_Z2(move(matrix_for_elimination));
    store_X_ijk();
    determine_sign_ambiguities();
  }
  
  int rank_over_Z2() const {
    return independent_rows_over_Z2.size();
  }
  int rank_over_Q() const {
    return independent_rows_over_Q.size();
  }

  bool are_derivations_traceless() const {
    return derivations_traceless;
  }
  vector<WeightAndCoefficient>&& weights_and_coefficients() && {
      return move(weights);
  }
 	list<SignConfiguration>&& sign_configurations() && {
  	return move(sign_configurations_);
  }
  Matrix M_Delta() const {
    return Matrix{weight_matrix};   
  }
  exvector X_ijk() const {
    return X_ijk_;
  }
  exvector kernel_of_MDelta_transpose() const {
    return solve_over_Q(complete_matrix_transposed_M_Delta(0),generate_variables<StructureConstant>(N.x,weight_matrix.rows()));
  }
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
		matrix nikolayevsky_derivation() const {
			int m=weight_matrix.rows(), n=weight_matrix.cols();			
			auto b=weight_matrix.mul(weight_matrix.transpose()).solve(matrix_of_unknowns(m),matrix_of_ones(m));
			auto v=weight_matrix.transpose()*b+matrix_of_ones(n);
			return ex_to<matrix>(v.evalm());
		}
};

#endif
