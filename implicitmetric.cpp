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
#include "implicitmetric.h"
#include "weightmatrix.h"
#include "linearinequalities.h"


template<typename NumType>
auto equals_zero () {
	static auto lambda= [](NumType x) {return x==0;};
	return lambda;
}

ImplicitMetric::ImplicitMetric(const string& name, const exvector& X_ijk) : 
	name_{name}, X_ijk{X_ijk},	in_coordinate_hyperplane{any_of(X_ijk.begin(),X_ijk.end(),equals_zero<ex>())} 
	{}




exvector integer_multiple(exvector v) {
	return v;
}



ex raise_each_and_multiply (const exvector& v,const exvector& coefficients) {
//	return transform_reduce(v.begin(),v.end(),1,coefficients.begin(),[](ex x, ex n) {return pow(abs(x),n);},std::multiplies<ex>{}).normal();
	ex product=1;
	auto j=coefficients.begin();
	for (auto i=v.begin();i!=v.end();++i,++j)
		product*=pow(*i,*j);
	return product.normal();
}

exvector polynomial_equations_for_existence (const exvector& X, const exvector& csquared) {
	auto basis_of_KerMDeltaTranspose=basis_from_generic_element<Unknown>(X);
	exvector equations;
	transform(basis_of_KerMDeltaTranspose.begin(),basis_of_KerMDeltaTranspose.end(),back_inserter(equations),
		[&X,&csquared] (auto& coeff) {return abs(raise_each_and_multiply(X,coeff))==
			raise_each_and_multiply(csquared,coeff)
		;});			
	return equations;
}

exvector ImplicitMetric::polynomial_equations_for_existence(const exvector& csquared) const {
	if (in_coordinate_hyperplane) return {}; else return ::polynomial_equations_for_existence(X_ijk,csquared);
}


template<typename Unknown>
bool affinespace_intersects_orthant(vector<Z2> epsilon, exvector generic_vector) {
	auto mul = [] (Z2 x, ex y) {return x.to_Z_star()*y;};
	transform(epsilon.begin(),epsilon.end(),generic_vector.begin(), generic_vector.begin(),mul);
	return LinearInequalities{generic_vector.begin(),generic_vector.end(),Unknown{}}.has_solution();
}


//condition 2
list<pair<SignConfiguration,SignConfiguration>> signatures_compatible_with_X(const WeightMatrix& weight_matrix,const exvector& X)
{
	assert(weight_matrix.rows()==0 || !X.empty());
	list<pair<SignConfiguration,SignConfiguration>> result;
	for (auto epsilon : SignConfiguration::all_configurations(weight_matrix.cols())) {
		auto MDelta2epsilon=weight_matrix.image_of(sign_configuration_to_vector(weight_matrix.cols(),epsilon));
//		if (all_of(MDelta2epsilon.begin(),MDelta2epsilon.end(),[](Z2 z) {return z==0;})) continue;	//do not consider "trivial" sign configurations
			if (affinespace_intersects_orthant<Unknown>((MDelta2epsilon), X))
				result.push_back(make_pair(epsilon,vector_to_sign_configuration(MDelta2epsilon)));
	}
	return result;
}


/*
//condition 2''
list<SignConfiguration> signatures_compatible_with_X_and_kerMDelta_transpose(const WeightMatrix& weight_matrix,const exvector& X)
{
	assert(weight_matrix.rows()==0 || !X.empty());
	list<SignConfiguration> result;
	auto MDeltatranspose=to_matrix(transpose(weight_matrix.M_Delta()));
	for (auto epsilon : SignConfiguration::all_configurations(weight_matrix.cols())) {
		auto MDelta2epsilon=weight_matrix.image_of(sign_configuration_to_vector(weight_matrix.cols(),epsilon));
//		if (all_of(MDelta2epsilon.begin(),MDelta2epsilon.end(),[](Z2 z) {return z==0;})) continue;	//do not consider "trivial" sign configurations
			if (affinespace_intersects_orthant<Unknown>((MDelta2epsilon), X) && sign_change_leaves_ker_invariant(MDelta2epsilon, MDeltatranspose))			
				result.push_back(epsilon);
	}
	return result;
}*/



DiagonalMetric::DiagonalMetric(const string& name,const WeightMatrix& weight_matrix, const exvector& X_ijk) : ImplicitMetric{name,X_ijk},	signatures{signatures_compatible_with_X(weight_matrix,X_ijk)} {}	


ostream& operator<<(ostream& os,SignConfiguration sign_configuration) {
	string sep;
	for (int i=0;i<sign_configuration.size();++i) {
		os<<sep;
		sep=" ";
		if (sign_configuration[i]>0) os<<"+"; else os<<"-";
	}
	return os;
}
