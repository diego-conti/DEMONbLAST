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


DiagramProperties:: DiagramProperties(const WeightMatrix& weight_matrix) : 
  rank_over_Z2{weight_matrix.rank_over_Z2()},
  X_ijk_{X_solving_nonRicciflat_Einstein(weight_matrix)},  
  nilsoliton_Y_ijk{X_solving_nilsoliton(weight_matrix)}
{
	nikolayevsky=to_matrix(transpose(weight_matrix.M_Delta())).image_of(nilsoliton_Y_ijk);
	for (auto& x: nikolayevsky) ++x;
	if (!LinearInequalities(nilsoliton_Y_ijk.begin(),nilsoliton_Y_ijk.end(),Unknown{}).has_solution()) nilsoliton_Y_ijk.clear();
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
  return sstream.str();
}


string signature(const exvector& metric) {
  int positive=count_if(metric.begin(),metric.end(),[](ex x) {return x>0;});
  int negative=metric.size()-positive;
  return "("+to_string(positive)+","+to_string(negative)+")";
}

string DiagramProperties::diagram_data() const {
 		stringstream sstream;
	  if (are_all_derivations_traceless()) {
	    sstream<<"derivations are traceless"<<endl;
	    sstream<<"X="<<horizontal(X_ijk_)<<endl;
   		if (is_X_ijk_in_coordinate_hyperplane()) sstream<<"X_ijk in coordinate hyperplane; "<<endl;
   	  else sstream<<"X_ijk not in coordinate hyperplane; "<<endl;
    }
    else {
	    sstream<<"Nikolayevsky derivation: "<<nikolayevsky<<endl;
	    if (!nilsoliton_Y_ijk.empty()) sstream<<"nilsoliton vector: "<<nilsoliton_Y_ijk<<endl;
	    else sstream<<"nilsoliton vector not >0"<<endl;
	  }
    sstream<<"rank over Z_2 = "<<rank_over_Z2<<endl;
    return sstream.str();
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

