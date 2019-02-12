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
#include "labeled_tree.h"
#include "partitions.h"
#include "weightbasis.h"
#include "parsetree.h"

list<LabeledTree> LabeledTree::labelings() const  {
		auto node=node_with_less_unlabeled_incoming_arrows();
		if (node==NO_NODE) return {*this};
		else return labelings_with_hint(node);
}

list<LabeledTree> LabeledTree::labelings_with_hint(int hint_node) const {
		list<LabeledTree> result;
		auto unlabeled_arrows = unlabeled_arrows_pointing_to(hint_node);
		auto arrow=unlabeled_arrows.begin();
		auto first=*arrow;
		while (++arrow!=unlabeled_arrows.end()) {
			auto labelings=labelings_with_bracket(hint_node, first, *arrow);
			result.splice(result.end(),labelings);	
		}
		return result;
}

LabeledTree LabeledTree::with_added_label(int to_node, int from_1, int from_2) const {
   auto arrows_with_added_label = arrows();	
   auto add_label = [to_node,from_1,from_2](LabeledArrow& arrow) {
		if (arrow.node_in==from_1 && arrow.node_out==to_node) arrow.label=from_2;
		else if (arrow.node_in==from_2 && arrow.node_out==to_node) arrow.label=from_1;    
   };
   for_each(arrows_with_added_label.begin(), arrows_with_added_label.end(), add_label);
	return {*this, move(arrows_with_added_label)};
}

list<LabeledTree> LabeledTree::labelings_with_bracket(int to_node, int from_1, int from_2) const {
	//check that from_1 does not already have an outgoing arrow with label from_2; equivalently, that the bracket [from_1,from_2] does not already appear
	auto is_same_bracket=[from_1,from_2] (const LabeledArrow& arrow) {
 		return arrow.node_in==from_1 && arrow.label==from_2;
 	};
   if (any_of(arrows().begin(),arrows().end(),is_same_bracket)) return {};
   else return with_added_label(to_node,from_1,from_2).labelings(); 
}

vector<int> LabeledTree::unlabeled_arrows_pointing_to(int node) const {
	vector<int> result;
	for (auto arrow: arrows()) 
		if (arrow.node_out==node && !arrow.has_label()) result.push_back(arrow.node_in);	
	return result;
}
	
vector<int> LabeledTree::number_of_unlabeled_incoming_arrows() const {
	vector<int> unlabeled_incoming_arrows(no_nodes);
	for (auto arrow: arrows()) if (!arrow.has_label()) ++unlabeled_incoming_arrows[arrow.node_out];		
	return unlabeled_incoming_arrows;
}
  
const WeightBasisAndProperties& LabeledTree::weight_basis(DiagramDataOptions diagram_data_options) const {
  if (weight_basis_==nullptr || !weight_basis_->matches(diagram_data_options)) weight_basis_=make_unique<WeightBasisAndProperties>(*this,diagram_data_options);
  return *weight_basis_;
}

LabeledTree::~LabeledTree()=default;

LabeledTree::LabeledTree(const Tree& tree) : TreeBase<LabeledArrow> (tree, tree.arrows()) {}
 
LabeledTree::LabeledTree(const LabeledTree& labeled_tree) : TreeBase<LabeledArrow>{labeled_tree} {}	
 
template<typename Iterator> Iterator position_of_minimal_nonzero_element(Iterator begin, Iterator end) {
	typename Iterator::value_type min_val=0;
	Iterator min_pos=end;
	for (Iterator i=begin;i!=end;++i) 
		if (*i!=0)
			if (min_val==0 || *i<min_val) {min_pos=i; min_val=*i;}
	return min_pos;
}

int LabeledTree::node_with_less_unlabeled_incoming_arrows() const {
	if (number_of_nodes()==0) return NO_NODE;
   auto unlabeled_incoming_arrows=number_of_unlabeled_incoming_arrows();
	auto minimal_element=position_of_minimal_nonzero_element(unlabeled_incoming_arrows.begin(),unlabeled_incoming_arrows.end());
	return minimal_element!=unlabeled_incoming_arrows.end()? minimal_element-unlabeled_incoming_arrows.begin() : NO_NODE;
}

void LabeledTree::canonicalize_order_increasing() {
	static auto compare = [](LabeledArrow arrow1, LabeledArrow arrow2) {
		return arrow1.node_out<arrow2.node_out || arrow1.node_out==arrow2.node_out &&	min(arrow1.node_in,arrow1.label) < min(arrow2.node_in,arrow2.label);
	};
	sort_arrows(compare);
}
void LabeledTree::canonicalize_order_decreasing() {
	static auto compare = [](LabeledArrow arrow1, LabeledArrow arrow2) {
		return arrow1.node_out>arrow2.node_out || arrow1.node_out==arrow2.node_out &&	max(arrow1.node_in,arrow1.label) > max(arrow2.node_in,arrow2.label);
	};
	sort_arrows(compare);
}

LabeledTree LabeledTree::from_stream(istream& s) {
	using namespace DiagramParser;
	int nodes=parse_number(s);
	LabeledTree tree{nodes,nodes};
	LabeledArrow arrow;
	while (parse_arrow(s,arrow)) {
		tree.add_arrow(arrow);	
		swap(arrow.node_in,arrow.label);
		tree.add_arrow(arrow);	
	}
	return tree;
}

string LabeledTree::as_string() const {
	stringstream s;
	s<<number_of_nodes();
	char sep=':';
	for (auto weight: weights()) {
		s<<sep<<' '<<weight.node_in1+1<<"->["<<weight.node_in2+1<<"]"<<weight.node_out+1;
		sep=',';
	}
	return s.str();
}


//3-combination of integers with repetitions
class Triplet {
	friend ostream& operator<<(ostream& os, Triplet triplet);
	int n1,n2,n3;
public:
	Triplet(int i,int j, int k) : n1{i}, n2{j}, n3{k} {
		if (n1>n2) swap(n1,n2);
		if (n1>n3) swap(n1,n3);
		if (n2>n3) swap(n2,n3);
	}
	bool operator==(Triplet other) const {
		return n1==other.n1 && n2==other.n2 && n3==other.n3;
	}
	bool operator<(Triplet other) const {
		if (n1<other.n1) return true;
		else if (n1==other.n1) {
			if (n2<other.n2) return true;
			else if (n2==other.n2 && n3<other.n3) return true;
		}
		return false;	
	}
	bool has_repetition() const {  
		return n1==n2 || n2==n3;
	}
};

ostream& operator<<(ostream& os, Triplet triplet) {
	return os<<"{"<<triplet.n1+1<<","<<triplet.n2+1<<","<<triplet.n3+1<<"}";
}


class DoubleIncomingArrows {
  friend ostream& operator<<(ostream&, const DoubleIncomingArrows&);
	using DoubleArrowsWithMultiplicity = map<Triplet,int>;
	vector<DoubleArrowsWithMultiplicity> double_incoming_arrows;
	void add_concatenated_arrow(Triplet double_arrow, int out) {
		if (!double_arrow.has_repetition())			//ignore double arrows of the form j->^i k->^i h
			++(double_incoming_arrows[out][double_arrow]);
	}
	void add_concatenated_arrow(Weight i, Weight j) {
		if (i.node_out == j.node_in1) add_concatenated_arrow({i.node_in1, i.node_in2, j.node_in2},j.node_out);
		else if (i.node_out == j.node_in2) add_concatenated_arrow({i.node_in1, i.node_in2, j.node_in1},j.node_out);
	}
public:
	DoubleIncomingArrows(const LabeledTree& tree) : double_incoming_arrows{tree.number_of_nodes()} {
		auto weights = tree.weights();
		for (auto i=weights.begin();i!=weights.end();++i) {
			auto j=i;
			while (++j!=weights.end()) {
				add_concatenated_arrow(*i,*j);	
				add_concatenated_arrow(*j,*i);	
			}
		}
	}
	bool has_double_incoming_arrow_with_multiplicity_three() {
   	  for (auto arrow: double_incoming_arrows)
	    for (auto pair : arrow)
	      if (pair.second>2) {
	        nice_log<<"double incoming arrow "<<pair.first<<" has multiplicity "<<pair.second<<endl; 
	        return true;
	        }
      return false;
	}
	
	bool formal_jacobi() const {		
		static auto has_multiplicity_greater_than_one = [](pair<Triplet,int> double_arrow_with_multiplicity) {return double_arrow_with_multiplicity.second>1;};
		static auto all_multiplicities_greater_than_one=
			 [](const DoubleArrowsWithMultiplicity& arrows) {return all_of(arrows.begin(),arrows.end(),has_multiplicity_greater_than_one);};
   return all_of(double_incoming_arrows.begin(),double_incoming_arrows.end(),all_multiplicities_greater_than_one);
	}
};

ostream& operator<<(ostream& os, const DoubleIncomingArrows& d) {
		for (int i=0;i<d.double_incoming_arrows.size();++i) {
			os<<"node "<<i+1<<endl;
			for (auto double_arrow_with_multiplicity : d.double_incoming_arrows[i])
				os<<double_arrow_with_multiplicity.first<<" * "<<double_arrow_with_multiplicity.second<<endl;
		}
	 return os;
}



bool satisfies_formal_jacobi(const LabeledTree& tree) {
	return DoubleIncomingArrows{tree}.formal_jacobi();	
}
bool has_double_incoming_arrow_with_multiplicity_three(const LabeledTree& tree) {
  return DoubleIncomingArrows{tree}.has_double_incoming_arrow_with_multiplicity_three();
}


