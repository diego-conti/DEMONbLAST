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
#include "adinvariantobstruction.h"


template<typename Iterator, typename Iterator2>
list<set<typename Iterator::value_type>> subsets_of(Iterator begin, Iterator2 end) {
	list<set<typename Iterator::value_type>> result;
	if (begin==end) return {{}};
	else {
		auto first=*begin;
		list<set<typename Iterator::value_type>> without_first=subsets_of(++begin,end);
		list<set<typename Iterator::value_type>> with_first;
		for (auto x: without_first) {
			x.insert(first);
			with_first.push_back(x);
		}
		without_first.splice(without_first.end(),with_first);
		return without_first;
	}
}

//a subspace of a nice Lie algebra spanned by nice basis elements
class NiceSubspace {
	set<int> generators;
public:
	template<typename T>	NiceSubspace(T&& generators) : generators{std::forward<T>(generators)} {}
	NiceSubspace(int n) : generators{consecutive_numbers(0,n)} {}
	NiceSubspace() =default;	
	static list<NiceSubspace> all_subspaces(const SetOfNodes& diagram) {
		vector<int> all_nodes(diagram.number_of_nodes());
		iota(all_nodes.begin(),all_nodes.end(),0);
		list<NiceSubspace> result;
		for (auto generators :  subsets_of(all_nodes.begin(),all_nodes.end())) result.push_back(generators);
		return result;
	}
	
	int dimension() const {return generators.size();}

	NiceSubspace reached_from(const LabeledTree& diagram) const {
		NiceSubspace result;
	  for (auto& arrow: diagram.arrows()) 
      if (generators.find(arrow.node_in)!=generators.end()) result.generators.insert(arrow.node_out);
		return result;
	}

	NiceSubspace commuting_with(const LabeledTree& diagram) const {
		NiceSubspace result{diagram.number_of_nodes()};	
	  for (auto& arrow: diagram.arrows()) 
      if (generators.find(arrow.node_in)!=generators.end()) {result.generators.erase(arrow.label);}
		return result;
	}
	
	NiceSubspace ending_in(const LabeledTree& diagram) const {
		NiceSubspace result{diagram.number_of_nodes()};	
	  for (auto& arrow: diagram.arrows()) 
      if (generators.find(arrow.node_out)==generators.end()) {result.generators.erase(arrow.node_in);}
    return result;
   }
   
   string to_string() const {return horizontal(generators," ");}
};

vector<int> generalized_lcs(const LabeledTree& diagram, NiceSubspace V) {
	vector<int> result;
	result.push_back(V.dimension());	
	while (V.dimension()>0) {
		V=V.reached_from(diagram);
		result.push_back(V.dimension());	
	}	
	return result;
}

vector<int> generalized_ucs(const LabeledTree& diagram, NiceSubspace V) {
	vector<int> result;
	NiceSubspace C_kV=V.commuting_with(diagram);
	result.push_back(0);
	result.push_back(C_kV.dimension());
	while (C_kV.dimension()!=diagram.number_of_nodes()) {
		C_kV=C_kV.ending_in(diagram);
		result.push_back(C_kV.dimension());	
	}
	return result;
}

int element_at_or_last_element_in(const vector<int>& v, int index) {
	return index<v.size()? v[index] : v.back();
}

bool passes_obstruction_for_ad_invariant_metric(const LabeledTree& diagram, const NiceSubspace& V) {
	auto lcs=	generalized_lcs(diagram,V);
	auto ucs= generalized_ucs(diagram,V);	
	int i=1;
	while (i<lcs.size() || i<ucs.size()) {
		auto l=	element_at_or_last_element_in(lcs,i);
		auto u=	element_at_or_last_element_in(ucs,i);		
		if (l+u!=diagram.number_of_nodes()) {
//				cout<<"Diagram "<<diagram.as_string()<<", obsturction failed for "<<V.to_string()<<endl;		
//				cout<<lcs<<endl<<ucs<<endl;
			return false;
		}
		++i;
	}
//	cout<<"Diagram "<<diagram.as_string()<<", obsturction passed for "<<V.to_string()<<endl;
	return true;
}


bool passes_obstruction_for_ad_invariant_metric(const LabeledTree& diagram) {
	auto all=NiceSubspace::all_subspaces(diagram);
	return all_of(all.begin(),all.end(),[&diagram] (const NiceSubspace& V) {return passes_obstruction_for_ad_invariant_metric(diagram,V);});
}

