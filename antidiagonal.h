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
#ifndef ANTIDIAGONAL_H
#define ANTIDIAGONAL_H

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


struct Couple {
	int first;
	int second;
	operator bool() const {return first!=second;}
	Couple& advance(int n) {
		if (second<n-1) ++second;
		else if (first<n-2) {
			++first;
			second=first+1;
		}
		else first=second=0;
		return *this;
	}
	bool overlaps(Couple x) {
		return first==x.first || second==x.first||first==x.second || second==x.second;
	}
};

class OrderTwoAutomorphism {
	list<Couple> couples;
	int n;
	
	bool compatible_from(list<Couple>::iterator i) {
		for (;i!=couples.end();++i) {
			auto j=i;
			while (++j!=couples.end())
				if (i->overlaps(*j)) return false;
		}
		return true;
	}
	
	bool advance_till_compatible_to_the_right(list<Couple>::iterator i) {
		do {
			i->advance(n);
		}
		while (*i && !compatible_from(i));
		return *i;	
	}
	bool advance_till_compatible(list<Couple>::iterator i) {
		if (!advance_till_compatible_to_the_right(i)) return false;
		for (auto j=i;j!=couples.begin();) {
			*--j=*i;
			if (!advance_till_compatible_to_the_right(j)) return false;
		}
		return true;
	}
	static int apply(int k, Couple trasposition) {
		if (k==trasposition.first) k= trasposition.second;
		else if (k==trasposition.second) k= trasposition.first;
		return k;
	}
public:
	OrderTwoAutomorphism() : n{0} {};
	OrderTwoAutomorphism(int n, int k) : n{n} {
		assert(n>=2*k);
		while (k--) {
			couples.push_front({0,1});
			if (!compatible_from(couples.begin())) advance_till_compatible(couples.begin());
		}		
	}
	OrderTwoAutomorphism(int n,std::initializer_list<Couple> list) : couples{list},n{n} {}
	OrderTwoAutomorphism& operator++() {
		auto iter=couples.begin();	
		assert(iter!=couples.end());
		while (!advance_till_compatible(iter)) {
			if (++iter==couples.end()) {
				couples.clear();
				break;
			}
		}
		return *this;
	}
	operator bool() const {return !couples.empty();}
	Weight apply(Weight weight) const {
		for (auto& couple : couples) {
			weight.node_in1=apply(weight.node_in1,couple);
			weight.node_in2=apply(weight.node_in2,couple);
			weight.node_out=apply(weight.node_out,couple);
		}
		return weight;
	}
	int apply(int node) const {
		for (auto& couple : couples) node=apply(node,couple);
		return node;	
	}
	string to_string() const {
		string result;
		for (auto c: couples) result+="("+std::to_string(c.first+1)+" "+std::to_string(c.second+1)+")";
		if (couples.empty()) result="()";
		return result;	
	}
	matrix to_matrix() const {
		matrix sigma(n,n);
		for (int i=0;i<n;++i) sigma(i,i)=1;
		for (auto couple: couples) {
			sigma(couple.first,couple.second)= sigma(couple.second,couple.first)=1;
			sigma(couple.first,couple.first)= sigma(couple.second,couple.second)=0;
		}
		return sigma;
	}
	matrix sigma_diagonal_metric(const exvector& coefficients) const {
		matrix sigma(n,n);
		for (int i=0;i<n;++i) sigma(i,i)=coefficients[i];
		for (auto couple: couples) {
			sigma(couple.first,couple.second)= sigma(couple.second,couple.first)=coefficients[couple.first];
			sigma(couple.first,couple.first)= sigma(couple.second,couple.second)=0;
		}
		return sigma;
	}	
	
};

inline ostream& operator<<(ostream& os, const OrderTwoAutomorphism& sigma) {
	return os<<sigma.to_string();
}

list<OrderTwoAutomorphism> ricci_flat_sigma(const WeightMatrix& weight_matrix);
bool sigma_defines_ricci_flat(const WeightMatrix& weight_matrix,const OrderTwoAutomorphism& sigma);

bool sigma_enhanced_has_cycle_of_length_two(const WeightMatrix& weight_matrix,const OrderTwoAutomorphism& sigma);

#endif
