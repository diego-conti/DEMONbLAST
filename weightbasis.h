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

class SignConfiguration {
	vector<int> signs;
public:
	SignConfiguration(int sign_ambiguities) {
		signs.insert(signs.begin(),sign_ambiguities,1);
	}
	bool operator==(const SignConfiguration& other) const {return signs==other.signs;}
	bool operator!=(const SignConfiguration& other) const {return signs!=other.signs;}
	SignConfiguration& operator++()  {
		int i=0;
			do {
				signs[i]=-signs[i];
		} while (signs[i]==1 && ++i<signs.size()) ;
		return *this;
	}
	bool has_next() const {
		return any_of(signs.begin(),signs.end(),[](int sign) {return sign==1;});
	}
	static list<SignConfiguration> all_configurations(int signs) {
		list<SignConfiguration> result;
		SignConfiguration conf{signs}; 
		result.push_back(conf); 
		while (conf.has_next())	result.push_back(++conf);
		return result;
	}
	int size() const {return signs.size();}
	int& operator[] (int i) {return signs[i];}	
	int operator[] (int i) const {return signs[i];}	
	template<typename T>
	auto multiply(vector<T> coefficients) const {
		assert(coefficients.size()>=signs.size());
		for (int i=0;i<signs.size();++i)
			coefficients[i]*=signs[i];
		return coefficients;
	}
	auto begin() const {return signs.begin();}
	auto end() const {return signs.end();}
};

ostream& operator<<(ostream& os,SignConfiguration sign_configuration);


ostream& operator<<(ostream& os,WeightAndCoefficient weight);

class WeightMatrix;

class DiagramProperties {
  int rank_over_Z2;
  bool X_ijk_in_coordinate_hyperplane;
  exvector X_ijk_;
  exvector nikolayevsky;
	exvector nilsoliton_Y_ijk;
	list<SignConfiguration> nilsoliton_signatures_;
	vector<WeightAndCoefficient> weights;	//this is only used to keep track of the order of the weights
public:
  DiagramProperties(const WeightMatrix& weight_matrix);
  bool are_all_derivations_traceless() const {return !X_ijk_.empty();}
  bool is_X_ijk_in_coordinate_hyperplane() const {return X_ijk_in_coordinate_hyperplane;}
  const exvector& X_ijk() const {return X_ijk_;}
 	virtual bool is_M_Delta_surjective() const=0;
 	virtual string diagram_data() const;
 	const	list<SignConfiguration>& signatures() const {return nilsoliton_signatures_;}
};

class DiagramPropertiesNonSurjectiveMDelta : public DiagramProperties {
  int no_rows;
  int rank_over_Q;
  bool X_ijk_in_coordinate_hyperplane;
  exvector kernel_of_MDelta_transpose;
public:
  string diagram_data() const ;
  DiagramPropertiesNonSurjectiveMDelta(const WeightMatrix& weight_matrix);
  bool is_M_Delta_surjective() const override {return false;}
};

class DiagramPropertiesSurjectiveMDelta : public DiagramProperties {
  list<exvector> einstein_metrics_;
  int rank_over_Q;

  list<exvector> compute_einstein_metrics(const Matrix& complete_matrix) const;  
public:
  DiagramPropertiesSurjectiveMDelta(const WeightMatrix& weight_matrix);
  bool is_M_Delta_surjective() const override {return true;}
  const list<exvector>& einstein_metrics() const {return einstein_metrics_;}
  bool is_X_ijk_in_reachable_orthant() const {
    return !(einstein_metrics_.empty());
   }
	string diagram_data() const override;
};



class WeightBasis {
	vector<WeightAndCoefficient> weights;
	list<SignConfiguration> configurations;	
  unique_ptr<DiagramProperties> diagram_properties;
  int number_of_nodes_;	
public:
	explicit WeightBasis(const LabeledTree& tree); 
	vector<WeightAndCoefficient> weights_and_coefficients() const {return weights;}
	list<SignConfiguration> sign_configurations() const {return configurations;}	
	const DiagramProperties& properties() const {return *diagram_properties;}
	int number_of_nodes() const {return number_of_nodes_;}
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

class EinsteinCoefficientConfiguration : public CoefficientConfiguration {
public:
	EinsteinCoefficientConfiguration(const WeightBasis& weight_basis) : CoefficientConfiguration{weight_basis.sign_configurations(),weight_basis.number_of_nodes()} {	
		assert(weight_basis.properties().are_all_derivations_traceless());
		auto& X_ijk=weight_basis.properties().X_ijk();
		auto coeff=X_ijk.begin();
		for (auto weight: weight_basis.weights_and_coefficients()) 
		    add_weight(weight,sqrt(abs(*coeff++)));  
  }
};
#endif
