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
#ifndef IMPLICIT_METRIC_H
#define IMPLICIT_METRIC_H

#include "includes.h"
#include <ginac/ginac.h>
#include "xginac.h"
#include <algorithm>
#include "gauss.h"
#include <wedge/liegroup.h>
#include "horizontal.h"
#include "log.h"


using namespace GiNaC;
using namespace Wedge;
using std::optional;
using std::nullopt;

class WeightMatrix;

class SignConfiguration {
	vector<int> signs;
public:
	SignConfiguration(int sign_ambiguities) {
		signs.insert(signs.begin(),sign_ambiguities,1);
	}
	template<typename T>
	SignConfiguration(const vector<T>& signs) : signs{signs} {}
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
	static list<SignConfiguration> all_configurations(int signs,int limit) {		
		list<SignConfiguration> result;
		SignConfiguration conf{signs}; 
		result.push_back(conf); 
 		if ((1<<signs)>limit) 
				nice_log<<"SignConfiguration::all_configurations invoked with signs="<<signs<<"; ignoring all sign configurations after the first "<<limit<<endl;
		while (conf.has_next() && result.size()<limit)	result.push_back(++conf);
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

class ImplicitMetric {
	string name_;
	exvector X_ijk;
	bool in_coordinate_hyperplane;
protected:
	ImplicitMetric()=default;
	ImplicitMetric(const string& name, const exvector& X_ijk);
	virtual void dump_extra(ostream& os) const {}
public:
	string to_string() const {
		stringstream sstream;
		sstream<<name_<<endl;
		sstream<<"X: "<<horizontal(X_ijk)<<endl;
		if (in_coordinate_hyperplane) sstream<<"in coordinate hyperplane"<<endl;
	  else dump_extra(sstream);
		return sstream.str();
	}
	const exvector& X() const {return X_ijk;}
	exvector polynomial_equations_for_existence(const exvector& csquared) const;	//return empty if there is no metric, the Lie algebra is abelian or M_Delta is surjective
	virtual string classification(const exvector& csquared) const {return "cannot classify";}
	virtual string solution_to_polynomial_equations_or_empty_string(const exvector& csquared) const {return {};}	
	bool is_in_coordinate_hyperplane() const {return in_coordinate_hyperplane;}
	string name() const {return name_;}
	virtual bool no_metric_regardless_of_polynomial_conditions() const {return is_in_coordinate_hyperplane();}
};

class DiagonalMetric : public ImplicitMetric {
	list<pair<SignConfiguration,SignConfiguration>> potential_signatures;	//first element of each pair is the signature, the other is M_Delta of it.
	optional<SignConfiguration> sign_configuration_from_image(const SignConfiguration& image) const;
	optional<set<vector<int>>> exact_signatures_for_codimension_one(const exvector& csquared) const;
	pair<bool,set<vector<int>>> riemannian_like_signatures() const;
	optional<set<vector<int>>> exact_signatures(const exvector& csquared) const;

	int dimension_coker_MDelta;
protected:
	void dump_extra(ostream& os) const override {
			if (potential_signatures.empty()) os<<"no metric (any signature)"<<endl;
			for (auto x: potential_signatures) 
			  os<<"(potential) signature: ("<<x.first<<") -> ("<<x.second<<")"<<endl;
	}
public:
	DiagonalMetric(const string& name, const WeightMatrix& weight_matrix, const exvector& X_ijk);
	bool no_metric_regardless_of_polynomial_conditions() const override {return is_in_coordinate_hyperplane()|| potential_signatures.empty();}
	string solution_to_polynomial_equations_or_empty_string(const exvector& csquared) const override;	
	string classification(const exvector& csquared) const override;
};

class SigmaCompatibleMetric : public ImplicitMetric {
	vector<int> automorphism;
protected:
	void dump_extra(ostream& os) const override {
		os<<"sigma="<<horizontal(automorphism)<<endl;
	}
public:	
	SigmaCompatibleMetric()=default;
	SigmaCompatibleMetric(const string& name, const WeightMatrix& weight_matrix, const exvector& X_ijk,const vector<int>& automorphism) : ImplicitMetric{name + ":"+horizontal(incremented(automorphism)), X_ijk}, automorphism{automorphism} {}
};

#endif
