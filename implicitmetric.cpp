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


template<typename Unknown> 
class PolynomialSolver {
	exvector solutions;	
	bool can_solve_;
	void solve_no_variables(ex equation) {
		can_solve_=true;
		if (equation.is_zero())
			solutions.push_back(lst{});
	}
	void solve_in_variable(ex equation, ex x) {
		auto a=equation.coeff(x,2);
		auto b=equation.coeff(x,1);
		auto c=equation.coeff(x,0);
		auto delta=b*b-4*a*c;
		switch (equation.degree(x)) {
			case 1: 
				can_solve_=true;
				solutions.push_back(x==-c/b);
				break;
			case 2: 
				can_solve_=true;
				if (delta.is_zero()) solutions.push_back(x==-b/(2*a));
				else if (delta>0) {
					solutions.push_back(x==(-b+sqrt(delta))/(2*a));
					solutions.push_back(x==(-b-sqrt(delta))/(2*a));
				}
				break;
			default:
					can_solve_=false;
			}
	}	
	PolynomialSolver()=default;
public:	
	PolynomialSolver(ex equation) {
		list<ex> variables;
		GetSymbols<Unknown> (variables,equation);
		if (variables.empty()) solve_no_variables(equation);
		else if (variables.size()>1) can_solve_=false;
		else solve_in_variable(equation,variables.front());
	}
	static PolynomialSolver solver_that_cannot_solve() {
		PolynomialSolver result;
		result.can_solve_=false;
		return result;	
	}
	bool can_solve() const {return can_solve_;}
	bool has_solution() const {return solutions.size();}
	list<ex> substitute_solutions(ex expression) const {
		list<ex> result;
		transform(solutions.begin(),solutions.end(),back_inserter(result),[expression] (ex solution) {return expression.subs(solution);});
		return result;
	}
	list<exvector> substitute_solutions(const exvector& expression) const {
		auto subs_into_vector=[&expression] (ex solution) {
			exvector result;
			transform(expression.begin(),expression.end(),back_inserter(result),[solution] (ex x) {return x.subs(solution);});
			return result;
		};	
		list<exvector> result;
		transform(solutions.begin(),solutions.end(),back_inserter(result),subs_into_vector);
		return result;
	}

	void add_solutions_from(const PolynomialSolver& solver) {
		if (!solver.can_solve()) can_solve_=false;
		solutions.insert(solutions.end(),solver.solutions.begin(),solver.solutions.end());
	}	
};



SignConfiguration sign_configuration_from_vector(const exvector& X) {
	vector<int> signs;
	transform(X.begin(),X.end(),back_inserter(signs),[](ex x) {return x>0? 1 : -1;});
	return {signs};
}


PolynomialSolver<Unknown> solutions_to_polynomial_problem_for_codimension_one(const exvector&  X,const exvector& csquared) {
	auto basis_of_KerMDeltaTranspose=basis_from_generic_element<Unknown>(X);
	assert (basis_of_KerMDeltaTranspose.size()==1);
	auto coeff=basis_of_KerMDeltaTranspose.front();
	ex lhs=raise_each_and_multiply(X,coeff);
	ex rhs=raise_each_and_multiply(csquared,coeff);
	auto solutions=PolynomialSolver<Unknown>{(lhs-rhs).numer()};
	solutions.add_solutions_from(PolynomialSolver<Unknown>{(lhs+rhs).numer()});
	return solutions;
}

vector<int> negative_signs(const SignConfiguration& sign_configuration) {
	vector<int> result;
	for (int i=0;i<sign_configuration.size();++i)
		if (sign_configuration[i]<0) result.push_back(i);
	return result;
}

string negative_signs_to_string(const vector<int>& negative_signs) {
	if (negative_signs.empty()) return "\\emptyset";
	else return horizontal(incremented(negative_signs),"");
}

string sign_configurations_to_string(const set<vector<int>>& negative_signs) {
	stringstream s;
	s<<"S=\\{";
		if (!negative_signs.empty()) {
			auto i=negative_signs.begin();
			s<<negative_signs_to_string(*i);
			while (++i!=negative_signs.end())
						s<<","<<negative_signs_to_string(*i);
		}
	s<<"\\}";
	return s.str();
}

optional<SignConfiguration> DiagonalMetric::sign_configuration_from_image(const SignConfiguration& image) const {
	for (auto& pair : potential_signatures)
		if (pair.second==image) return pair.first;
	return nullopt;
}

optional<set<vector<int>>> DiagonalMetric::exact_signatures_for_codimension_one(const exvector& csquared) const {
	set<vector<int>> signatures;
	auto solutions=solutions_to_polynomial_problem_for_codimension_one(X(),csquared);
	if (!solutions.can_solve()) return nullopt;
	for (auto Xsol : solutions.substitute_solutions(X())) {
		auto conf=sign_configuration_from_vector(Xsol);
		auto signature=sign_configuration_from_image(conf);
		if (signature) signatures.insert(negative_signs(signature.value()));
	}
	return signatures;
}

string DiagonalMetric::classification(const exvector& csquared) const {
	optional<set<vector<int>>> signatures;
	if (no_metric_regardless_of_polynomial_conditions())
			signatures.emplace();	//assign empty set: no valid signature
	else if (dimension_coker_MDelta==0) {
			signatures.emplace();
			for (auto& signature : potential_signatures) signatures.value().insert(negative_signs(signature.first));
	}
	else if (dimension_coker_MDelta==1)  
			signatures=exact_signatures_for_codimension_one(csquared);
	return  signatures? sign_configurations_to_string(signatures.value()) : ImplicitMetric::classification(csquared);
}

//TODO fai scrivere le segnature
string DiagonalMetric::solution_to_polynomial_equations_or_empty_string(const exvector& csquared) const {
	stringstream s;
	if (dimension_coker_MDelta==0) {
		s<<"X="<<horizontal(X());
		s<<"-> ("<<sign_configuration_from_vector(X())<<")"<<endl;
	}
	else if (dimension_coker_MDelta==1)	{
		auto solutions=solutions_to_polynomial_problem_for_codimension_one(X(),csquared);
		if (solutions.can_solve()) {
			for (auto Xsol : solutions.substitute_solutions(X())) {
				s<<"X="<<horizontal(Xsol);
				s<<"-> ("<<sign_configuration_from_vector(Xsol)<<")"<<endl;
			}		
			if (!solutions.has_solution()) s<<"no solution"<<endl;
		}
	}
	return s.str();
}


DiagonalMetric::DiagonalMetric(const string& name,const WeightMatrix& weight_matrix, const exvector& X_ijk) :
	ImplicitMetric{name,X_ijk},	
	potential_signatures{signatures_compatible_with_X(weight_matrix,X_ijk)},
	dimension_coker_MDelta{weight_matrix.rows()-weight_matrix.rank_over_Q()}
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
