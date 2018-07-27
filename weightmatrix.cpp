#include "weightmatrix.h"


template<typename Builder>
Matrix complete_with_constant_vector(Builder&& builder, ex lambda) {
 	return adjoin(builder,ConstantMatrixBuilder{builder.rows(),1,lambda});  
}

exvector X_solving_nonRicciflat_Einstein(const WeightMatrix& weight_matrix) {	
	return solve_over_Q(complete_with_constant_vector(transpose(weight_matrix.M_Delta()),1),generate_variables<StructureConstant>(N.x,weight_matrix.rows()));
};
exvector X_solving_Ricciflat(const WeightMatrix& weight_matrix) {	
	return solve_over_Q(complete_with_constant_vector(transpose(weight_matrix.M_Delta()),0),generate_variables<StructureConstant>(N.x,weight_matrix.rows()));
};

exvector X_solving_nilsoliton(const WeightMatrix& weight_matrix) {	
	auto MDelta=weight_matrix.M_Delta();
	auto complete_matrix = adjoin(transpose(MDelta),ConstantMatrixBuilder{weight_matrix.cols(),1,1},ColumnVectorMatrixBuilder{nikolayevsky_derivation(MDelta)});
	exvector sol= solve_over_Q(move(complete_matrix),generate_variables<Unknown>(N.x,weight_matrix.rows()+1));
	sol.pop_back();		
	return sol;
};

exvector nikolayevsky_derivation(const Matrix& MDelta) {
	auto gram = complete_with_constant_vector(matrix_product(MDelta,transpose(MDelta)),1);
	auto b=	solve_over_Q(gram,generate_variables<Unknown>(N.x,MDelta.rows()));
	auto v=to_matrix(transpose(MDelta)).image_of(b);
	for (auto& x: v) ++x;
	return v;
}

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

WeightMatrix::WeightMatrix(vector<WeightAndCoefficient> unordered_weights,int dimension) : weights(unordered_weights.size()), cols_{dimension} { 
	int last=unordered_weights.size();
	for (auto row : ::independent_rows_over_Z2(matrix_from_weights(unordered_weights,dimension))) {
			unordered_weights[row].eliminate_parameter();
			weights[--last]=unordered_weights[row];
	 		++independent_rows_over_Z2;
	}
	for (auto parameter : ReorderedMatrix{unordered_weights,dimension}.independent_rows_over_Q()) {
	 	++independent_rows_over_Q;
	  if (!unordered_weights[parameter].parameter_eliminated()) {
	 		unordered_weights[parameter].eliminate_parameter();
			weights[--last]=unordered_weights[parameter];
		}		
	}
	for (auto& weight: unordered_weights) 
		if (!weight.parameter_eliminated())
			weights[--last]=weight;
	assert(last==0);
}
