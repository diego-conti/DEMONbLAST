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
using std::multiset;

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


exvector multiply_signs (const vector<Z2>& epsilon, exvector X) {
	auto mul = [] (Z2 x, ex y) {return x.to_Z_star()*y;};
	transform(epsilon.begin(),epsilon.end(),X.begin(), X.begin(),mul);
	return X;
}

//condition 2'' (unused!)
list<pair<SignConfiguration,SignConfiguration>> signatures_compatible_with_X_and_kerMDelta_transpose(const WeightMatrix& weight_matrix,const exvector& X)
{
	assert(weight_matrix.rows()==0 || !X.empty());
	auto MDelta=weight_matrix.M_Delta();
	auto gram = to_matrix( matrix_product(MDelta,transpose(MDelta)));
	exvector ones(MDelta.rows(),ex{1});
	list<pair<SignConfiguration,SignConfiguration>> result;
	for (auto epsilon : SignConfiguration::all_configurations(weight_matrix.cols())) {
		auto MDelta2epsilon=weight_matrix.image_of(sign_configuration_to_vector(weight_matrix.cols(),epsilon));
		auto epsilonX=multiply_signs(MDelta2epsilon,X);
//		if (all_of(MDelta2epsilon.begin(),MDelta2epsilon.end(),[](Z2 z) {return z==0;})) continue;	//do not consider "trivial" sign configurations
		if (affinespace_intersects_orthant<Unknown>((MDelta2epsilon), X) && gram.image_of(epsilonX)==ones)
			result.push_back(make_pair(epsilon,vector_to_sign_configuration(MDelta2epsilon)));
	}
	return result;
}


/*
//condition 2'' (some old version)
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

//TODO write  a function that attempts to solve a generic polynomial equation (and fails if degree exceeds 2)
list<ex> solve_second_degree_equation(ex eq) {
	list<ex> variables;
	GetSymbols<Unknown> (variables,eq);
	if (variables.empty()) {
		ex x=Unknown{N.x};
		if (eq.is_zero()) return {x==x}; else return {};
	}
	else if (variables.size()>1) {
			std::cerr<<eq<<endl;
			throw std::invalid_argument("equation with more than one variable");
	}
	auto x=variables.front();
	auto a=eq.coeff(x,2);
	auto b=eq.coeff(x,1);
	auto c=eq.coeff(x,0);
	auto delta=b*b-4*a*c;
	switch (eq.degree(x)) {
		case 1: return {x==-c/b};
		case 2: 
			if (delta.is_zero()) return {x==-b/(2*a)} ;
			else if (delta<0) return {};
			else return  {x==(-b+sqrt(delta))/(2*a),x==(-b-sqrt(delta))/(2*a)};
		default:
				std::cerr<<eq<<endl;
			 throw std::invalid_argument("not a polynomial of degree 1 or 2");
		}
}

SignConfiguration sign_configuration_from_vector(const exvector& X) {
	vector<int> signs;
	transform(X.begin(),X.end(),back_inserter(signs),[](ex x) {return x>0? 1 : -1;});
	return {signs};
}

string DiagonalMetric::solution_to_polynomial_equations_or_empty_string(const exvector& csquared) const {
	stringstream s;
	auto basis_of_KerMDeltaTranspose=basis_from_generic_element<Unknown>(X());
	if (basis_of_KerMDeltaTranspose.size()!=1) return {};
	auto coeff=basis_of_KerMDeltaTranspose.front();
	multiset<ex,ex_is_less> as_multiset{coeff.begin(),coeff.end()};
	if (as_multiset.count(1)==2 && as_multiset.count(-1)==2 && as_multiset.count(0)==coeff.size()-4) {
		ex lhs=raise_each_and_multiply(X(),coeff);
		ex rhs=raise_each_and_multiply(csquared,coeff);
		try {
			auto solutions=solve_second_degree_equation((lhs-rhs).numer());
			solutions.splice(solutions.end(),	solve_second_degree_equation((lhs+rhs).numer()));
			if (solutions.empty()) s<<"no solution"<<endl;
			else for (auto sol: solutions) {
				exvector Xsol;
				transform(X().begin(),X().end(),back_inserter(Xsol),[&sol] (ex x) {return x.subs(sol);});
				s<<"X="<<horizontal(Xsol);
				s<<"-> ("<<sign_configuration_from_vector(Xsol)<<")"<<endl;
			}
		}
		catch (std::exception& e) {
			std::cerr<<e.what()<<endl;
			std::cerr<<X()<<endl;
			std::cerr<<csquared<<endl;
			std::cerr<<coeff<<endl;

			std::cerr<<lhs<<endl;
			std::cerr<<rhs<<endl;
			std::cerr<<(lhs-rhs).numer()<<endl;
			std::cerr<<(lhs+rhs).numer()<<endl;
			throw;
		}
	}
	return s.str();
}


DiagonalMetric::DiagonalMetric(const string& name,const WeightMatrix& weight_matrix, const exvector& X_ijk) :
	ImplicitMetric{name,X_ijk},	
	potential_signatures{signatures_compatible_with_X(weight_matrix,X_ijk)}
	 {}	


ostream& operator<<(ostream& os,SignConfiguration sign_configuration) {
	string sep;
	for (int i=0;i<sign_configuration.size();++i) {
		os<<sep;
		sep=" ";
		if (sign_configuration[i]>0) os<<"+"; else os<<"-";
	}
	return os;
}
