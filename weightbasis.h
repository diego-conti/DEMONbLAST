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

#include <iostream>
#include <list>
#include <ginac/ginac.h>
#include "xginac.h"
#include <algorithm>
#include "gauss.h"
#include <wedge/liegroup.h>
#include "labeled_tree.h"
#include "options.h"
#include "implicitmetric.h"
#include "antidiagonal.h"
#include "diagramanalyzer.h"

using namespace std;
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
	bool parameter_eliminated() const {return parameter_eliminated_;}
	bool same_weight_as(Weight weight) const {
		return weight.node_in1==node_in1 && weight.node_in2==node_in2 && weight.node_out==node_out;
	}
private:
	bool parameter_eliminated_=false;
};

ostream& operator<<(ostream& os,WeightAndCoefficient weight);

class WeightMatrix;

class DiagramProperties {
public:
  DiagramProperties(const WeightMatrix& weight_matrix, const list<vector<int>>& automorphisms, DiagramDataOptions options);
 	virtual bool is_M_Delta_surjective() const=0;
 	virtual string diagram_data() const;
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
		}
 		return result; 		
 	}
  bool are_all_derivations_traceless() const {
		return all_of(nikolayevsky.begin(),nikolayevsky.end(),[] (ex x) {return x.is_zero();});
	}
	bool has_nontrivial_automorphisms() const {return !automorphisms.empty();}
	
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
protected:
	DiagramDataOptions options;
	virtual void print_matrix_data(ostream& os) const;
private:
	static string RICCI_FLAT_DIAGONAL() {return "Ricci-flat(diag)"s;}
	static string RICCI_FLAT_SIGMA() {return "Ricci-flat(sigma)"s;}
	static string NILSOLITON_DIAGONAL() {return "nilsoliton(diag)"s;}
  int rank_over_Z2;
  list<unique_ptr<ImplicitMetric>> metrics;
	vector<WeightAndCoefficient> weights;	//this is only used to keep track of the order of the weights TODO canonicalize the order and do away with this field
	exvector nikolayevsky;
	list<vector<int>> automorphisms;
	string imMDelta2;
	list<OrderTwoAutomorphism> ricci_flat_antidiagonal;
	DiagramAnalyzer diagram_analyzer;		//combinatorial data about the diagram
};

class DiagramPropertiesNonSurjectiveMDelta : public DiagramProperties {
  int no_rows;
  int rank_over_Q;
  bool X_ijk_in_coordinate_hyperplane;
  exvector kernel_of_MDelta_transpose;
protected:
 	void print_matrix_data(ostream& os) const override;	
public:
  DiagramPropertiesNonSurjectiveMDelta(const WeightMatrix& weight_matrix,const list<vector<int>>& automorphisms, DiagramDataOptions options);
  bool is_M_Delta_surjective() const override {return false;}
};

class DiagramPropertiesSurjectiveMDelta : public DiagramProperties {
  int rank_over_Q;
protected:
 	void print_matrix_data(ostream& os) const override;	
public:
  DiagramPropertiesSurjectiveMDelta(const WeightMatrix& weight_matrix, const list<vector<int>>& automorphisms, DiagramDataOptions options);
  bool is_M_Delta_surjective() const override {return true;}
	string diagram_data() const override;
};


class WeightBasis {
	list<SignConfiguration> configurations;	
  int number_of_nodes_;	
protected:
	vector<WeightAndCoefficient> weights;
	explicit WeightBasis(const WeightMatrix& weight_matrix, const list<vector<int>>& automorphisms); 
public:
	explicit WeightBasis(const LabeledTree& tree); 
	vector<WeightAndCoefficient> weights_and_coefficients() const {return weights;}
	list<SignConfiguration> sign_configurations() const {return configurations;}	
	int number_of_nodes() const {return number_of_nodes_;}
};

class WeightBasisAndProperties : public WeightBasis {
  unique_ptr<DiagramProperties> diagram_properties;
	explicit WeightBasisAndProperties(const WeightMatrix& weight_matrix, const list<vector<int>>& automorphisms, DiagramDataOptions diagram_data_options); 
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


struct WeightAndValue : Weight {
	WeightAndValue(Weight weight, ex coefficient) : Weight(weight), value{coefficient} {}
	ex value;
	WeightAndValue operator*=(int sign) {value*=sign; return *this;}
};

class CoefficientConfiguration {
public:
	CoefficientConfiguration(const CoefficientConfiguration&)=delete;
	int lie_algebra_dimension() const {return nodes;}
//iteration through all configurations
	CoefficientConfiguration& operator++() {
		++current;
		return *this;
	}	
	operator bool () const {
		return current!=sign_configurations.end();
	}	
//iterations through this configuration
	vector<WeightAndValue> weights() const {
		return current->multiply(configuration_with_positive_signs);
	}
protected:
  void add_weight(WeightAndCoefficient weight, ex coefficient) {
      configuration_with_positive_signs.emplace_back(weight,coefficient);
  }
  CoefficientConfiguration(const list<SignConfiguration> & sign_configurations, int nodes) : sign_configurations{sign_configurations}, nodes{nodes} {
  	if (sign_configurations.empty()) this->sign_configurations.emplace_back(0);
  	current=this->sign_configurations.begin();
  }
private:
	vector<WeightAndValue> configuration_with_positive_signs;
	list<SignConfiguration> sign_configurations;
	list<SignConfiguration>::const_iterator current;
	const int nodes;
};

class CoefficientConfigurationWithoutRedundantParameter : public CoefficientConfiguration {
public:
	CoefficientConfigurationWithoutRedundantParameter(const WeightBasis& weight_basis) : CoefficientConfiguration{weight_basis.sign_configurations(),weight_basis.number_of_nodes()} {
		int no_parameters=0;
		nice_log<<weight_basis.weights_and_coefficients()<<endl;
		for (auto weight: weight_basis.weights_and_coefficients()) 
			if (weight.parameter_eliminated()) add_weight(weight,1);
			else add_weight(weight,StructureConstant{N.a(++no_parameters)});
  }
};

enum class MetricType {
	NONFLAT_NILSOLITON,RICCIFLAT
};

class MetricCoefficientConfiguration : public CoefficientConfiguration {
	const exvector& X(const WeightBasisAndProperties& weight_basis, MetricType metric_type) {
		switch (metric_type) {
			case MetricType::NONFLAT_NILSOLITON:
				return weight_basis.properties().diagonal_nilsoliton_metric()->X();
			case MetricType::RICCIFLAT:
				return weight_basis.properties().diagonal_ricci_flat_metric()->X();
			default:
				throw 0;
			}
		}
public:
	MetricCoefficientConfiguration(const WeightBasisAndProperties& weight_basis, MetricType metric_type) : CoefficientConfiguration{weight_basis.sign_configurations(),weight_basis.number_of_nodes()} {	
		auto& X_ijk=X(weight_basis,metric_type);
		auto coeff=X_ijk.begin();
		for (auto weight: weight_basis.weights_and_coefficients()) 
		    add_weight(weight,sqrt(abs(*coeff++)));  
  }
};

#endif
