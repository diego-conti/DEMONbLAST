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
#ifndef COEFFICIENT_CONFIGURATION_H
#define COEFFICIENT_CONFIGURATION_H

#include "weightbasis.h"

struct WeightAndValue : Weight {
	WeightAndValue(Weight weight, ex coefficient) : Weight(weight), value{coefficient} {}
	ex value;
	WeightAndValue operator*=(int sign) {value*=sign; return *this;}
};

class CoefficientConfiguration {
public:
	CoefficientConfiguration(int nodes) : nodes{nodes} {}
	CoefficientConfiguration(const CoefficientConfiguration&)=delete;
	int lie_algebra_dimension() const {return nodes;}
//iteration through all configurations
	virtual CoefficientConfiguration& operator++()=0;
	virtual operator bool () const=0;
//iterations through this configuration. The resulting vector may not have the same order as the vector passed in the constructor
	virtual vector<WeightAndValue> weights() const =0;
	virtual ~CoefficientConfiguration()=default;
private:
	const int nodes;
};


class CoefficientConfigurationWithVariableSigns : public CoefficientConfiguration {
public:
//iteration through all configurations
	CoefficientConfigurationWithVariableSigns& operator++() override {
		++current;
		return *this;
	}	
	operator bool () const override {
		return current!=sign_configurations.end();
	}	
//iterations through this configuration. The resulting vector may not have the same order as the vector passed in the constructor
	vector<WeightAndValue> weights() const override {
		return current->multiply(configuration_with_positive_signs);
	}
protected:
	struct Weights {
		vector<WeightAndValue> weights_with_sign_eliminated;
		vector<WeightAndValue> weights_with_sign;
		vector<WeightAndValue> join() &&  {
			weights_with_sign.insert(weights_with_sign.end(),weights_with_sign_eliminated.begin(),weights_with_sign_eliminated.end());
			return move(weights_with_sign);
		}
		void add_weight(WeightAndCoefficient weight, ex coefficient) {
			if (weight.sign_and_parameter_eliminated()) 
	      weights_with_sign_eliminated.emplace_back(weight,coefficient);
  		else
	      weights_with_sign.emplace_back(weight,coefficient);
		}
	};
  
  CoefficientConfigurationWithVariableSigns(const list<SignConfiguration> & sign_configurations, int nodes,Weights&& weights) : CoefficientConfiguration{nodes},
	 	configuration_with_positive_signs{move(weights).join()}, sign_configurations{sign_configurations}
	{
  	if (sign_configurations.empty()) this->sign_configurations.emplace_back(0);
  	current=this->sign_configurations.begin();
  }
private:
	vector<WeightAndValue> configuration_with_positive_signs;
	list<SignConfiguration> sign_configurations;
	list<SignConfiguration>::const_iterator current;
};

class CoefficientConfigurationWithoutRedundantParameter : public CoefficientConfigurationWithVariableSigns {
	static Weights weights_and_values(const WeightBasis& weight_basis) {
		Weights weights;
		int no_parameters=0;
		nice_log<<weight_basis.weights_and_coefficients()<<endl;
		for (auto weight: weight_basis.weights_and_coefficients()) 
			if (weight.parameter_eliminated()) weights.add_weight(weight,1);
			else weights.add_weight(weight,StructureConstant{N.a(++no_parameters)});
		return weights;	
	}
public:
	CoefficientConfigurationWithoutRedundantParameter(const WeightBasis& weight_basis) : 
		CoefficientConfigurationWithVariableSigns{weight_basis.sign_configurations(),weight_basis.number_of_nodes(),weights_and_values(weight_basis)} {}
};

enum class MetricType {
	NONFLAT_NILSOLITON,RICCIFLAT
};

class MetricCoefficientConfiguration : public CoefficientConfigurationWithVariableSigns {
	static const exvector& X(const WeightBasisAndProperties& weight_basis, MetricType metric_type) {
		switch (metric_type) {
			case MetricType::NONFLAT_NILSOLITON:
				return weight_basis.properties().diagonal_nilsoliton_metric()->X();
			case MetricType::RICCIFLAT:
				return weight_basis.properties().diagonal_ricci_flat_metric()->X();
			default:
				throw 0;
			}
		}
	static Weights weights_and_values(const WeightBasis& weight_basis, const exvector& X) {
		Weights weights;
		auto coeff=X.begin();
		nice_log<<weight_basis.weights_and_coefficients()<<endl;
		for (auto weight: weight_basis.weights_and_coefficients()) 
			weights.add_weight(weight,sqrt(abs(*coeff++)));  
		return weights;	
	}

public:
	MetricCoefficientConfiguration(const WeightBasisAndProperties& weight_basis, MetricType metric_type) 
		: CoefficientConfigurationWithVariableSigns{weight_basis.sign_configurations(),weight_basis.number_of_nodes(),weights_and_values(weight_basis,X(weight_basis,metric_type))} {}
};

class CoefficientLists {
	vector<exvector> coefficient_lists; 
public:
	CoefficientLists()=default;
	CoefficientLists(const CoefficientLists& c)=default;
	CoefficientLists(std::initializer_list<exvector> c) : coefficient_lists{move(c)} {}
	auto begin() const {return coefficient_lists.begin();}
	auto end() const {return coefficient_lists.end();}	
};

class FixedCoefficientConfiguration : public CoefficientConfiguration {
	vector<Weight> weights_;
	CoefficientLists coefficients_lists;
	vector<exvector>::const_iterator current;
public:
	FixedCoefficientConfiguration(const WeightBasisAndProperties& weight_basis, const CoefficientLists& coefficients_lists)  : 
		CoefficientConfiguration{weight_basis.number_of_nodes()}, coefficients_lists{coefficients_lists} {
		for (auto& w: weight_basis.weights_and_coefficients()) weights_.push_back(w);
		current=this->coefficients_lists.begin();
	}

	FixedCoefficientConfiguration(const FixedCoefficientConfiguration&)=delete;
	FixedCoefficientConfiguration& operator++() override {
		++current;
		return *this;
	}	
	operator bool () const override {
		return current!=coefficients_lists.end();
	}	
	vector<WeightAndValue> weights() const override {
		vector<WeightAndValue> result;
		auto w=weights_.begin();
		for (auto c=current->begin();c!=current->end();++w,++c)		
			result.emplace_back(*w,*c);
		return result;
	}
};


#endif
