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
#include "weightbasis.h"
#include "liegroupsfromdiagram.h"
#include "horizontal.h"
#include "weightmatrix.h"
#include "linearinequalities.h"
#include "antidiagonal.h"

using namespace Wedge;
using namespace std;

ostream& operator<<(ostream& os,WeightAndCoefficient weight) {
	os<<"["<<weight.node_in1+1<<","<<weight.node_in2+1<<"]=";
	if (weight.parameter_eliminated()) 	os<<"";
	else os<<"\\pm";
	os<<weight.node_out+1;
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



bool has_order_two(const vector<int>& nontrivial_automorphism) {
	for (int i=0;i<nontrivial_automorphism.size();++i)
		if (nontrivial_automorphism[nontrivial_automorphism[i]]!=i) return false;
	return true;
}

exvector symmetrize(const exvector& X, const vector<int>& sigma) {
	assert(X.size()==sigma.size());
	exvector result(X.size());
	for (int i=0;i<X.size();++i)
		result[i]=X[i]+X[sigma[i]];
	return result;
}


unique_ptr<DiagonalMetric> make_diagonal_metric(const string& name, const WeightMatrix& weight_matrix, const exvector& X_ijk, DiagramDataOptions options) {
	return options.only_riemannian_like_metrics()? make_unique<DiagonalMetric>(name,weight_matrix,X_ijk,only_riemannian_like_metrics) :  make_unique<DiagonalMetric>(name,weight_matrix,X_ijk);
}

DiagramProperties:: DiagramProperties(const WeightMatrix& weight_matrix, const list<vector<int>>& automorphisms, DiagramDataOptions options) : 
 	options{options},
 	no_rows{weight_matrix.M_Delta().rows()},
 	no_cols{weight_matrix.M_Delta().cols()},
  rank_over_Q{weight_matrix.rank_over_Q()},
  rank_over_Z2{weight_matrix.rank_over_Z2()},
  automorphisms{automorphisms},
  imMDelta2{options.sign_configurations_limit? "not computed" : image_mod2(weight_matrix).to_string()}, 	
 	ricci_flat_antidiagonal{options.with_antidiagonal_ricci_flat_sigma()? ricci_flat_sigma(weight_matrix) : list<OrderTwoAutomorphism>{}},
 	diagram_analyzer{options.analyze_diagram()? DiagramAnalyzer{weight_matrix.cols(),vector<WeightAndCoefficient>{weight_matrix.weight_begin(),weight_matrix.weight_end()}} : DiagramAnalyzer{}},
  kernel_of_MDelta_transpose{X_solving_Ricciflat(weight_matrix)}
{
	auto nilsoliton_X=X_solving_nilsoliton(weight_matrix);
	b=nilsoliton_X;	
	B=accumulate(b.begin(),b.end(),ex{0});
	auto diagonal_ricci_flat_X=X_solving_Ricciflat(weight_matrix);
	if (options.with_diagonal_nilsoliton_metrics())
		metrics.push_back(make_diagonal_metric(NILSOLITON_DIAGONAL(), weight_matrix,nilsoliton_X,options));
	if (options.with_diagonal_ricci_flat_metrics())
		metrics.push_back(make_diagonal_metric(RICCI_FLAT_DIAGONAL(), weight_matrix,diagonal_ricci_flat_X,options));
	nikolayevsky=to_matrix(transpose(weight_matrix.M_Delta())).image_of(nilsoliton_X);
	for (auto& x: nikolayevsky) ++x;
	if (options.with_sigma_compatible_ricci_flat_metrics()) 
		for (auto& sigma: automorphisms)
			if (has_order_two(sigma))
				metrics.push_back(make_unique<SigmaCompatibleMetric> (RICCI_FLAT_SIGMA(),weight_matrix, symmetrize(diagonal_ricci_flat_X,weight_matrix.sigma_on_VDelta(sigma)), sigma));
}

  
string signature(const exvector& metric) {
  int positive=count_if(metric.begin(),metric.end(),[](ex x) {return x>0;});
  int negative=metric.size()-positive;
  return "("+to_string(positive)+","+to_string(negative)+")";
}




void DiagramProperties::print_matrix_data(ostream& os) const {
	  if (are_all_derivations_traceless()) {
	    os<<"derivations are traceless"<<endl;
    }
    else {
	    os<<"Nikolayevsky derivation: "<<nikolayevsky<<endl;
	  }
    os<<"rank over Z_2 = "<<rank_over_Z2<<endl;
	  os<<"rank over Q = "<<rank_over_Q<<endl;
	  os<<"dim ker M_Delta = "<< dimension_kernel_M_Delta()<<endl;
	  os<<"dim coker M_Delta = "<<dimension_cokernel_M_Delta()<<endl;
	  if (dimension_cokernel_M_Delta()) {
		  os<<" kernel of (M_Delta)^T "<<kernel_of_MDelta_transpose<<"; generators:"<<endl;
		 	for (auto v : basis_from_generic_element<Unknown>(kernel_of_MDelta_transpose)) os<<v<<endl;
		  list<ex> symbols;
		  for (auto entry : kernel_of_MDelta_transpose)
		  for (auto symbol : symbols) 
		  	if (abs(entry.coeff(symbol))>1) os<<"unexpected coefficient"<<endl;  
		 }
    os<<"B="<<B<<endl;
    os<<"b="<<horizontal(b)<<endl;
}


string DiagramProperties::diagram_data() const {
 		stringstream sstream;
 		if (options.with_matrix_data()) print_matrix_data(sstream);
 		if (options.with_im_delta2())
	    sstream<<"Im M_Delta2: "<<endl<<imMDelta2<<endl;
	  for (auto& metric : metrics) 
	   	sstream<<metric->to_string()<<endl;
	  if (options.analyze_diagram())
	  	sstream<<diagram_analyzer.analysis()<<endl;
    /*
    for (auto x: SignConfiguration::all_configurations(nilsoliton_Y_ijk.size())) {
			exvector with_signs;
			auto mul = [] (int x, ex y) {return x*y;};
			transform(x.begin(),x.end(),nilsoliton_Y_ijk.begin(), back_inserter(with_signs),mul);
			if ( LinearInequalities{with_signs.begin(),with_signs.end(),Unknown{}}.has_solution()) sstream<<"orthant : "<<x<<endl;
		}*/
		if (options.with_antidiagonal_ricci_flat_sigma()) {
	    if (ricci_flat_antidiagonal.empty()) sstream<<"no ricci-flat antidiagonal sigma"<<endl;
			else sstream<<"ricci-flat antidiagonal sigma: "<<cut_at(ricci_flat_antidiagonal,options.ricci_flat_antidiagonal_limit)<<endl;
		}
    return sstream.str();
}


WEDGE_DECLARE_NAMED_ALGEBRAIC(Parameter,realsymbol);


exvector subs(exvector v,ex subs) {
  for (auto& x: v) x=x.subs(subs);
  return v;
}

exvector minus_one_exp(exvector v) {
  for (auto& x: v)
      x=reduce_mod_Z2(x.expand()).is_zero()? 1: -1;
  return v;
}

/*
list<exvector> affine_space_over_Z2(const exvector& generic_element) {
  list<Parameter> parameters;
  GetSymbols<Parameter> (parameters,generic_element.begin(),generic_element.end());
  if (parameters.empty()) return {minus_one_exp(generic_element)};
  else {
    auto with_0=affine_space_over_Z2(subs(generic_element,parameters.front()==0));
    with_0.splice(with_0.end(),affine_space_over_Z2(subs(generic_element,parameters.front()==1)));
    return with_0;
  }
}*/

WeightBasis::WeightBasis(const WeightMatrix& weight_matrix,const list<vector<int>>& automorphisms, DiagramDataOptions diagram_data_options) : number_of_nodes_{weight_matrix.cols()} {
  configurations=diagram_data_options.sign_configurations_limit? 
  	SignConfigurations{weight_matrix,automorphisms,diagram_data_options.sign_configurations_limit}.sign_configurations() : SignConfigurations{weight_matrix,automorphisms}.sign_configurations();  
  weights.insert(weights.end(), weight_matrix.weight_begin(),weight_matrix.weight_end());
  if (!configurations.empty()) {
  		for (auto w: weights) nice_log<<w<<", ";
  		nice_log<<endl<<"remaining sign configurations:"<<endl;
			for (auto sc : configurations)
				nice_log<<sc<<endl;
	}
}
WeightBasis::WeightBasis(const LabeledTree& tree) : WeightBasis{WeightMatrix{tree.weights(),tree.number_of_nodes()},tree.nontrivial_automorphisms(),{}} {}

WeightBasisAndProperties::WeightBasisAndProperties(const WeightMatrix& weight_matrix, const list<vector<int>>& automorphisms, DiagramDataOptions options) 
	: WeightBasis{weight_matrix,automorphisms,options}
{
    diagram_properties = make_unique< DiagramProperties>(weight_matrix,automorphisms,options);
}

WeightBasisAndProperties::WeightBasisAndProperties(const WeightMatrix& weight_matrix,const LabeledTree& tree, DiagramDataOptions options)
  : WeightBasisAndProperties{weight_matrix,
      (options.with_automorphisms() || (options.use_automorphisms_to_eliminate_signs() && weight_matrix.rank_over_Z2()<weight_matrix.rows())) ? tree.nontrivial_automorphisms() : list<vector<int>>{},
          options
} {}

WeightBasisAndProperties::WeightBasisAndProperties(const LabeledTree& tree, DiagramDataOptions options) 
	: WeightBasisAndProperties{WeightMatrix{tree.weights(),tree.number_of_nodes()},tree,options} {}

EnhancedWeightBasis::EnhancedWeightBasis(const LabeledTree& tree, OrderTwoAutomorphism sigma) : WeightBasis{weight_matrix(tree,sigma), automorphisms(tree,sigma),{}} {}

WeightMatrix EnhancedWeightBasis::weight_matrix(const LabeledTree& tree, OrderTwoAutomorphism sigma) {
		list<Weight> weights=tree.weights();
		for (auto weight: weights) weights.push_front(sigma.apply(weight));
		return WeightMatrix{weights,tree.number_of_nodes()};
}
