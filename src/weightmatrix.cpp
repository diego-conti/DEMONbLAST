#include "weightmatrix.h"


template<typename Builder>
Matrix complete_with_constant_vector(Builder&& builder, ex lambda) {
 	return adjoin(builder,ConstantMatrixBuilder{builder.rows(),1,lambda});  
}

exvector X_solving_nonRicciflat_Einstein(const WeightMatrix& weight_matrix) {	
	return solve_over_Q(complete_with_constant_vector(transpose(weight_matrix.M_Delta()),1),generate_variables<Unknown>(N.x,weight_matrix.rows()));
};
exvector X_solving_Ricciflat(const WeightMatrix& weight_matrix) {	
	return solve_over_Q(complete_with_constant_vector(transpose(weight_matrix.M_Delta()),0),generate_variables<Unknown>(N.x,weight_matrix.rows()));
};

exvector X_solving_nilsoliton(const WeightMatrix& weight_matrix) {	
	auto MDelta=weight_matrix.M_Delta();
	auto gram = complete_with_constant_vector(matrix_product(MDelta,transpose(MDelta)),1);
	auto b=	solve_over_Q(gram,generate_variables<Unknown>(N.x,MDelta.rows()));
	return b;
};

struct ReorderedMatrix {
	Matrix m;
	vector<int> row_to_parameter_correspondence;
public:
	ReorderedMatrix(const vector<WeightAndCoefficient>& weights, int dimension) : m(weights.size(),dimension), row_to_parameter_correspondence(m.rows()) {	
 	  int next_row_in_front=0, next_row_in_back=m.rows();
		for (int i=0;i<weights.size();++i)
  		if (weights[i].parameter_eliminated()) {
  		  row_to_parameter_correspondence[next_row_in_front]=i;
  		  populate_row(next_row_in_front++,weights[i],m);
  	}
     else {
       	row_to_parameter_correspondence[--next_row_in_back]=i;
  		  populate_row(next_row_in_back,weights[i],m);
  	}
	}
	vector<int>	independent_rows_over_Q() && {
		auto reordered=::independent_rows_over_Q(move(m));
		vector<int> result(reordered.size());
		transform(reordered.begin(),reordered.end(),result.begin(),
			[this] (int row) {return row_to_parameter_correspondence[row];}
			);
		return result;
	}
};

Matrix matrix_from_weights(const vector<WeightAndCoefficient>& weights,int dimension) {
	Matrix result{weights.size(),dimension};
	for (int i=0;i<weights.size();++i) 
		populate_row(i,weights[i],result);
	return result;
}



WeightMatrix::WeightMatrix(vector<WeightAndCoefficient>&& unordered_weights,int dimension) : weights{move(unordered_weights)}, cols_{dimension} { 
	basis_over_Z2 = ::independent_rows_over_Z2(matrix_from_weights(weights,dimension));
	independent_rows_over_Z2=basis_over_Z2.size();
	for (auto row : basis_over_Z2) 
			weights[row].eliminate_sign_and_parameter();
	
	for (auto parameter : ReorderedMatrix{weights,dimension}.independent_rows_over_Q()) {
	 	++independent_rows_over_Q;
	  if (!weights[parameter].parameter_eliminated()) {
	 		weights[parameter].eliminate_parameter();
		}		
	}
	for (int i=0;i<weights.size();++i)
		if (!weights[i].sign_and_parameter_eliminated())
			complement_of_basis_over_Z2.push_back(i);
}

string sign_configuration_to_string(const SignConfiguration& sign_configuration) {
	stringstream s;
	s<<sign_configuration;
	return s.str();
}

vector<Z2> sign_configuration_to_vector(int dimension, const SignConfiguration& epsilon) {
	vector<Z2> w_epsilon(dimension);
	auto logsign = [] (int sign) {return sign>0? 0 : 1;};
	transform(epsilon.begin(),epsilon.end(),w_epsilon.begin(),logsign);
	return w_epsilon;
}

SignConfiguration vector_to_sign_configuration(const vector<Z2>& epsilon) {
	vector<int> as_integers(epsilon.size());
	transform(epsilon.begin(),epsilon.end(), as_integers.begin(), [] (Z2 x) {return x.to_Z_star();} );
	return as_integers;
}


ImageMod2 image_mod2(const WeightMatrix& weight_matrix) {
	ImageMod2 result;
	for (auto delta : SignConfiguration::all_configurations(weight_matrix.cols())) {
		auto vector=sign_configuration_to_vector(weight_matrix.cols(),delta);
		result.insert(sign_configuration_to_string(delta), weight_matrix.image_of(vector));
	}
	return result;
}

