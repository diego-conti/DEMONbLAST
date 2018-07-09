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
#ifndef WEIGHT_BASIS_H
#define WEIGHT_BASIS_H

#include "includes.h"
#include <ginac/ginac.h>
#include "xginac.h"
#include "gauss.h"
#include <wedge/liegroup.h>
#include "labeled_tree.h"
#include "options.h"
#include "implicitmetric.h"
#include "antidiagonal.h"
#include "diagramanalyzer.h"

using namespace GiNaC;
using namespace Wedge;


struct WeightAndCoefficient : Weight {
	WeightAndCoefficient() = default;
	WeightAndCoefficient(const WeightAndCoefficient&) = default;
	WeightAndCoefficient(WeightAndCoefficient&&) = default;
	WeightAndCoefficient(Weight&& weight) : Weight(weight) {};
	WeightAndCoefficient(const Weight& weight) : Weight(weight) {};	
	WeightAndCoefficient& operator=(const WeightAndCoefficient&)=default;
	WeightAndCoefficient& operator=(WeightAndCoefficient&&)=default;
	void eliminate_parameter() {parameter_eliminated_=true;}
	void eliminate_sign_and_parameter() {parameter_eliminated_=true; sign_eliminated_=true;}
	bool parameter_eliminated() const {return parameter_eliminated_;}
	bool sign_and_parameter_eliminated() const {return sign_eliminated_;}
	bool same_weight_as(Weight weight) const {
		return weight.node_in1==node_in1 && weight.node_in2==node_in2 && weight.node_out==node_out;
	}
private:
	bool parameter_eliminated_=false;
	bool sign_eliminated_=false;
};

ostream& operator<<(ostream& os,WeightAndCoefficient weight);

class WeightMatrix;

class DiagramProperties {
public:
  DiagramProperties(const WeightMatrix& weight_matrix, const list<vector<int>>& automorphisms, DiagramDataOptions options);
 	int dimension_kernel_M_Delta() const {return no_cols-rank_over_Q;}
 	int dimension_cokernel_M_Delta() const {return no_rows-rank_over_Q;}
 	string diagram_data() const;
 	bool potentially_admits_metrics() const {
 		assert(!metrics.empty());
 		auto result=any_of(metrics.begin(),metrics.end(),
 			[](auto& metric) {return !metric->no_metric_regardless_of_polynomial_conditions();});
 		return result;
 	}
 	string polynomial_conditions(const exvector& csquared) const {
 		string result;
 		for (auto& metric : metrics) {
 			if (metric->no_metric_regardless_of_polynomial_conditions()) continue;
 			result+=metric->name();
 			auto polynomial_equations=metric->polynomial_equations_for_existence(csquared);
			if (!polynomial_equations.empty())	result+=" when structure constants satisfy: "+horizontal(polynomial_equations)+"\n";
			else result+=" always\n";
			result +=metric->solution_to_polynomial_equations_or_empty_string(csquared);
		}
 		return result; 		
 	}
  bool are_all_derivations_traceless() const {
		return all_of(nikolayevsky.begin(),nikolayevsky.end(),[] (ex x) {return x.is_zero();});
	}
	bool has_nontrivial_automorphisms() const {return !automorphisms.empty();}
	
	string classification_of_metrics(const exvector& csquared) const {
		string result;
		for (auto& metric : metrics)
			result+=metric->name()+"&"+metric->classification(csquared);
		return result;
	}
	
	const ImplicitMetric* diagonal_nilsoliton_metric() const {
		auto it=find_if(metrics.begin(),metrics.end(),[] (auto& metric) {return metric->name()==NILSOLITON_DIAGONAL();});
		return it!=metrics.end()? it->get() : nullptr;
	}
	const ImplicitMetric* diagonal_ricci_flat_metric() const {
		auto it=find_if(metrics.begin(),metrics.end(),[] (auto& metric) {return metric->name()==RICCI_FLAT_DIAGONAL();});
		return it!=metrics.end()? it->get() : nullptr;
	}
	bool matches(DiagramDataOptions options) const {return this->options==options;}
	list<OrderTwoAutomorphism> automorphisms_giving_ricci_flat() const {return ricci_flat_antidiagonal;}
	exvector nikolayevsky_derivation() const {return nikolayevsky;}
	bool simple_nikolayevsky_derivation() const {return set<ex,ex_is_less>(nikolayevsky.begin(),nikolayevsky.end()).size()==nikolayevsky.size();}
private:
	static string RICCI_FLAT_DIAGONAL() {return "Ricci-flat(diag)"s;}
	static string RICCI_FLAT_SIGMA() {return "Ricci-flat(sigma)"s;}
	static string NILSOLITON_DIAGONAL() {return "nilsoliton(diag)"s;}
	DiagramDataOptions options;
	void print_matrix_data(ostream& os) const;
  int no_rows, no_cols, rank_over_Q;
  int rank_over_Z2;
  list<unique_ptr<ImplicitMetric>> metrics;
	exvector nikolayevsky;
	list<vector<int>> automorphisms;
	string imMDelta2;
	list<OrderTwoAutomorphism> ricci_flat_antidiagonal;
	DiagramAnalyzer diagram_analyzer;		//combinatorial data about the diagram
	ex B;
	exvector b;
  exvector kernel_of_MDelta_transpose;
};



class WeightBasis {
	list<SignConfiguration> configurations;	
  int number_of_nodes_;	
protected:
	vector<WeightAndCoefficient> weights;
	explicit WeightBasis(const WeightMatrix& weight_matrix, const list<vector<int>>& automorphisms, DiagramDataOptions diagram_data_options); 
public:
	explicit WeightBasis(const LabeledTree& tree); 
	vector<WeightAndCoefficient> weights_and_coefficients() const {return weights;}
	list<SignConfiguration> sign_configurations() const {return configurations;}	
	int number_of_nodes() const {return number_of_nodes_;}
};

class WeightBasisAndProperties : public WeightBasis {
  unique_ptr<DiagramProperties> diagram_properties;
	explicit WeightBasisAndProperties(const WeightMatrix& weight_matrix, const list<vector<int>>& automorphisms, DiagramDataOptions diagram_data_options); 
	explicit WeightBasisAndProperties(const WeightMatrix& weight_matrix,const LabeledTree& tree, DiagramDataOptions options);
public:
	explicit WeightBasisAndProperties(const LabeledTree& tree, DiagramDataOptions diagram_data_options); 
	const DiagramProperties& properties() const {return *diagram_properties;}
	bool matches(DiagramDataOptions options) const {return diagram_properties && diagram_properties->matches(options);}
};

class EnhancedWeightBasis  : public WeightBasis {
	WeightMatrix weight_matrix(const LabeledTree& tree, OrderTwoAutomorphism sigma);
	list<vector<int>> automorphisms(const LabeledTree& tree, OrderTwoAutomorphism sigma) {
		return tree.nontrivial_automorphisms(); //FIXME there will be more.
	}
public:
	EnhancedWeightBasis(const LabeledTree& tree, OrderTwoAutomorphism sigma);
};




template<typename Parameter>
vector<exvector> basis_from_generic_element(const exvector& generic_vector) {	
	list<ex> free; 
	GetSymbols<Parameter>(free, generic_vector.begin(),generic_vector.end());
	vector<exvector> basis;
	transform (free.begin(),free.end(),back_inserter(basis),
		[&generic_vector] (ex x) {return coefficients_of(generic_vector,x);}
	);
	return basis;
}

#endif
