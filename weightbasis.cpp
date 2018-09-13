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
#include "weightbasis.h"
#include <wedge/ginaclinalg.h>
#include <sstream>
#include "liegroupsfromdiagram.h"
#include "horizontal.h"
#include "weightmatrix.h"
#include "linearinequalities.h"

using namespace Wedge;


ostream& operator<<(ostream& os,WeightAndCoefficient weight) {
	os<<"["<<weight.node_in1+1<<","<<weight.node_in2+1<<"]=";
	if (weight.parameter_eliminated()) 	os<<"";
	else os<<"\\pm";
	os<<weight.node_out+1;
	return os;
}

ostream& operator<<(ostream& os,SignConfiguration sign_configuration) {
	for (int i=0;i<sign_configuration.size();++i)
		if (sign_configuration[i]>0) os<<"+ "; else os<<"- ";
	return os;
}


int logsign(ex x) {
  int result= (x>0)? 0 : 1;
  return result;
}
exvector logsign(const exvector& non_zero_numbers) {
  exvector result;
  for (ex x: non_zero_numbers) result.push_back(logsign(x));
  return result;
}

template<typename Unknown>
bool affinespace_intersects_orthant(vector<Z2> epsilon, exvector generic_vector) {
	auto mul = [] (Z2 x, ex y) {return x.to_Z_star()*y;};
	transform(epsilon.begin(),epsilon.end(),generic_vector.begin(), generic_vector.begin(),mul);
	return LinearInequalities{generic_vector.begin(),generic_vector.end(),Unknown{}}.has_solution();
}

bool sign_change_leaves_ker_invariant(const vector<Z2>& epsilon, const Matrix& matrix) {
	auto vector_in_kernel=	solve_over_Q(adjoin(matrix, ConstantMatrixBuilder{matrix.rows(),1,0}),generate_variables<Unknown>(N.x,matrix.cols()));
	auto mul = [] (Z2 x, ex y) {return x.to_Z_star()*y;};
	transform(epsilon.begin(),epsilon.end(),vector_in_kernel.begin(), vector_in_kernel.begin(),mul);
	auto image=matrix.image_of(vector_in_kernel);
	return all_of(image.begin(),image.end(),[](ex x) {return x.is_zero();});
}

Matrix wa_minus_p(const WeightMatrix& weight_matrix) {
	auto MDeltatranspose=to_matrix(transpose(weight_matrix.M_Delta()));
	auto b=X_solving_nilsoliton(weight_matrix);
	if (b.empty()) return MDeltatranspose;
	ex trace_b=accumulate(b.begin(),b.end(),ex{});	
	auto p=MDeltatranspose.image_of(b);
	assert (!trace_b.is_zero());
	for (int i=0;i<MDeltatranspose.rows();++i)
		for (int j=0;j<MDeltatranspose.cols();++j)
			MDeltatranspose(i,j)-=p[i]/trace_b;
	return MDeltatranspose;
}

//condition 2
list<SignConfiguration> signatures_compatible_with_X(const WeightMatrix& weight_matrix,const exvector& X)
{
	assert(weight_matrix.rows()==0 || !X.empty());
	list<SignConfiguration> result;
	for (auto epsilon : SignConfiguration::all_configurations(weight_matrix.cols())) {
		auto MDelta2epsilon=weight_matrix.image_of(sign_configuration_to_vector(weight_matrix.cols(),epsilon));
//		if (all_of(MDelta2epsilon.begin(),MDelta2epsilon.end(),[](Z2 z) {return z==0;})) continue;	//do not consider "trivial" sign configurations
			if (affinespace_intersects_orthant<Unknown>((MDelta2epsilon), X))
				result.push_back(epsilon);
	}
	return result;
}


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
}

exvector subs_in_vector(exvector v, ex subs) {
	transform(v.begin(),v.end(),v.begin(),[&subs] (ex x) {return x.subs(subs);});
	return v;
}

template<typename Parameter>
vector<exvector> basis_from_generic_element(const exvector& generic_vector) {	
	list<ex> free; 
	GetSymbols<Parameter>(free, generic_vector.begin(),generic_vector.end());
	list<ex> substitutions;
	transform(free.begin(),free.end(),back_inserter(substitutions),[](ex symbol) {return symbol==0;});
	vector<exvector> basis;
	transform (free.begin(),free.end(),back_inserter(basis),
		[&substitutions,&generic_vector] (ex x) {return subs_in_vector(subs_in_vector(generic_vector,x==1),lst{substitutions});}
	);
	return basis;
}

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

exvector polynomial_equations_for_existence (const exvector& ricci_flat_X, const exvector& csquared) {
	auto basis_of_KerMDeltaTranspose=basis_from_generic_element<Unknown>(ricci_flat_X);
	exvector equations;
	transform(basis_of_KerMDeltaTranspose.begin(),basis_of_KerMDeltaTranspose.end(),back_inserter(equations),
		[&ricci_flat_X,&csquared] (auto& coeff) {return abs(raise_each_and_multiply(ricci_flat_X,coeff))==
			raise_each_and_multiply(csquared,coeff)
		;});	
	return equations;
}


DiagramProperties:: DiagramProperties(const WeightMatrix& weight_matrix) : 
  rank_over_Z2{weight_matrix.rank_over_Z2()},
  X_ijk_{X_solving_nonRicciflat_Einstein(weight_matrix)},  
  nilsoliton_Y_ijk{X_solving_nilsoliton(weight_matrix)},
  ricci_flat_Y_ijk{X_solving_Ricciflat(weight_matrix)},
  nilsoliton_signatures_{signatures_compatible_with_X(weight_matrix, nilsoliton_Y_ijk)},
  ricci_flat_signatures_{signatures_compatible_with_X(weight_matrix, ricci_flat_Y_ijk)},
  weights{weight_matrix.weight_begin(),weight_matrix.weight_end()}
{
	nikolayevsky=to_matrix(transpose(weight_matrix.M_Delta())).image_of(nilsoliton_Y_ijk);
	for (auto& x: nikolayevsky) ++x;
  X_ijk_in_coordinate_hyperplane=any_of(begin(X_ijk_),end(X_ijk_),[](ex x) {return x.expand().is_zero();});
}

DiagramPropertiesNonSurjectiveMDelta:: DiagramPropertiesNonSurjectiveMDelta(const WeightMatrix& weight_matrix) 
  : DiagramProperties(weight_matrix), no_rows(weight_matrix.M_Delta().rows()),  rank_over_Q{weight_matrix.rank_over_Q()},
  kernel_of_MDelta_transpose(X_solving_Ricciflat(weight_matrix))
{}
  
string DiagramPropertiesNonSurjectiveMDelta::diagram_data() const {
  stringstream sstream;
	sstream<<DiagramProperties::diagram_data(); 		
  sstream<<"rank over Q = "<<rank_over_Q;
  sstream<< "< "<< no_rows<< "(M_Delta not surjective);"<<endl;
  sstream<<" kernel of (M_Delta)^T "<<kernel_of_MDelta_transpose<<endl;
  list<ex> symbols;
  GetSymbols<Unknown>(symbols,kernel_of_MDelta_transpose.begin(),kernel_of_MDelta_transpose.end());
  for (auto entry : kernel_of_MDelta_transpose)
  for (auto symbol : symbols) 
  	if (abs(entry.coeff(symbol))>1) sstream<<"unexpected coefficient"<<endl;  
  return sstream.str();
}


string signature(const exvector& metric) {
  int positive=count_if(metric.begin(),metric.end(),[](ex x) {return x>0;});
  int negative=metric.size()-positive;
  return "("+to_string(positive)+","+to_string(negative)+")";
}

string DiagramProperties::diagram_data() const {
 		stringstream sstream;
 		sstream<<horizontal(weights)<<endl;
	  if (are_all_derivations_traceless()) {
	    sstream<<"derivations are traceless"<<endl;
	    sstream<<"X="<<horizontal(X_ijk_)<<endl;
   		if (is_X_ijk_in_coordinate_hyperplane()) sstream<<"X_ijk in coordinate hyperplane; "<<endl;
   	  else sstream<<"X_ijk not in coordinate hyperplane; "<<endl;
    }
    else {
	    sstream<<"Nikolayevsky derivation: "<<nikolayevsky<<endl;
	  }
	  sstream<<"Nilsoliton vector [k=1]: "<<horizontal(nilsoliton_Y_ijk)<<endl;
	  if (nilsoliton_signatures_.empty()) sstream<<"no diagonal nilsoliton metric (any signature)"<<endl;
	  for (auto x: nilsoliton_signatures_) 
	  	sstream<<"(potential) nilsoliton signature: "<<x<<endl;	
	  if (ricci_flat_signatures_.empty()) sstream<<"no diagonal Ricci-flat metric (any signature)"<<endl;
  	for (auto x: ricci_flat_signatures_) 
	  	sstream<<"(potential) Ricci-flat signature: "<<x<<endl;	
    
    for (auto x: SignConfiguration::all_configurations(nilsoliton_Y_ijk.size())) {
			exvector with_signs;
			auto mul = [] (int x, ex y) {return x*y;};
			transform(x.begin(),x.end(),nilsoliton_Y_ijk.begin(), back_inserter(with_signs),mul);
			if ( LinearInequalities{with_signs.begin(),with_signs.end(),Unknown{}}.has_solution()) sstream<<"orthant : "<<x<<endl;
		}

    sstream<<"rank over Z_2 = "<<rank_over_Z2<<endl;
    return sstream.str();
}

exvector DiagramProperties::polynomial_equations_for_existence_of_ricci_flat_metric(const exvector& csquared) const {
  if (ricci_flat_signatures_.empty()) return {};
  else return polynomial_equations_for_existence(ricci_flat_Y_ijk,csquared);
}

string DiagramPropertiesSurjectiveMDelta::diagram_data() const  {
  stringstream sstream;
 		sstream<<DiagramProperties::diagram_data(); 		
	  sstream<<"rank over Q = "<<rank_over_Q<< "(M_Delta surjective);"<<endl; 		
 		if (is_X_ijk_in_reachable_orthant()) sstream<<"X_ijk in reachable orthant; "<<endl;
 	  else sstream<<"X_ijk not in reachable orthant; "<<endl;
 	  for (auto& einstein_metric : einstein_metrics()) sstream<<"Einstein metric: "<<einstein_metric<<", signature "<<signature(einstein_metric)<<endl;
    return sstream.str();
}

WEDGE_DECLARE_NAMED_ALGEBRAIC(Parameter,realsymbol);

DiagramPropertiesSurjectiveMDelta::DiagramPropertiesSurjectiveMDelta(const WeightMatrix& weight_matrix) 
  : DiagramProperties{weight_matrix},rank_over_Q{weight_matrix.rank_over_Q()}
{
  if (are_all_derivations_traceless() && !is_X_ijk_in_coordinate_hyperplane()) 
   einstein_metrics_ = compute_einstein_metrics(complete_matrix(weight_matrix.M_Delta(), logsign(X_ijk())));
}

exvector subs(exvector v,ex subs) {
  for (auto& x: v) x=x.subs(subs);
  return v;
}

exvector minus_one_exp(exvector v) {
  for (auto& x: v)
      x=reduce_mod_Z2(x.expand()).is_zero()? 1: -1;
  return v;
}


list<exvector> affine_space_over_Z2(const exvector& generic_element) {
  list<Parameter> parameters;
  GetSymbols<Parameter> (parameters,generic_element.begin(),generic_element.end());
  if (parameters.empty()) return {minus_one_exp(generic_element)};
  else {
    auto with_0=affine_space_over_Z2(subs(generic_element,parameters.front()==0));
    with_0.splice(with_0.end(),affine_space_over_Z2(subs(generic_element,parameters.front()==1)));
    return with_0;
  }
}


list<exvector> DiagramPropertiesSurjectiveMDelta::compute_einstein_metrics(const Matrix& complete_matrix) const {
    nice_log<<"DiagramPropertiesSurjectiveMDelta::compute_einstein_metric"<<endl;
    nice_log<<complete_matrix<<endl;
    exvector vars = generate_variables<Parameter>(N.x,complete_matrix.cols()-1);
    auto generic_einstein_metric = solve_over_Z2(complete_matrix,vars);   
    if (generic_einstein_metric.empty()) return {};
    return affine_space_over_Z2(generic_einstein_metric);     
}



WeightBasis::WeightBasis(const LabeledTree& tree) : number_of_nodes_{tree.number_of_nodes()} {
  WeightMatrix weight_matrix{tree.weights(),tree.number_of_nodes()};
  configurations=SignConfigurations{weight_matrix, tree.nontrivial_automorphisms()}.sign_configurations();
  weights.insert(weights.end(), weight_matrix.weight_begin(),weight_matrix.weight_end());
  if (!configurations.empty()) {
  		nice_log<<"remaining sign configurations:"<<endl;
			for (auto sc : configurations)
				nice_log<<sc<<endl;
	}
  if (weight_matrix.rank_over_Q()==weights.size())
    diagram_properties = make_unique< DiagramPropertiesSurjectiveMDelta>(weight_matrix);
  else 
    diagram_properties = make_unique< DiagramPropertiesNonSurjectiveMDelta>(weight_matrix);
}

