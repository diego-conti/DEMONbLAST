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
#include <algorithm>
#include <iostream>
#include <ginac/ginac.h>
#include "tree.h"
#include "permutations.h"
#include "partitions.h"
#include "labeled_tree.h"

using namespace GiNaC;

ostream& operator<<(ostream& os, Arrow arrow) {
  return 	os<< arrow.node_in+1 << " -> "<< arrow.node_out+1<<endl;
}
ostream& operator<<(ostream& os, LabeledArrow arrow) {
		os<< arrow.node_in+1 << " -> "<< arrow.node_out+1;
		if (arrow.has_label()) os<< "[label="<< arrow.label+1 << "]";
		return os<<endl;
}

ostream& operator<<(ostream& os, LabeledDoubleArrow arrow) {
		os<< arrow.node_in+1 << " -> "<< arrow.node_out+1;
		if (arrow.has_label()) os<< "[label=\""<< arrow.label.first+1<<","<<arrow.label.second+1<<"\"]";
		return os<<endl;
}

template<typename BooleanType>
bool negate_in_place(BooleanType& x) {return (x=!x);}

ListOfArrows& ListOfArrows::operator++() {
    auto i=arrow_present.begin();
    while (!negate_in_place(*i++));
    assert(i!=arrow_present.end() || is_this_the_end());
    return *this;
}

ListOfArrows::ListOfArrowsConstIterator& ListOfArrows::ListOfArrowsConstIterator::operator++() {
  while (!list_of_arrows.at(++index) && index<list_of_arrows.size());
  return *this;
}



template<typename Closure>
list<Tree> add_one_or_more_nodes(const Tree& tree, int nodes, ListOfArrows first_list_of_arrows, Closure&& is_valid_tree);

template<typename Closure>
list<Tree> add_nodes(const Tree& tree, int nodes,  const ListOfArrows&  first_list_of_arrows,  Closure&& is_valid_tree) {
	if (nodes==0) return is_valid_tree(tree)? list<Tree>{tree} : list<Tree>{} ;
  else return add_one_or_more_nodes(tree,nodes, first_list_of_arrows,std::forward<Closure>(is_valid_tree));
}

template<typename Closure>
list<Tree> add_one_or_more_nodes(const Tree& tree, int nodes,  ListOfArrows first_list_of_arrows,  Closure&& is_valid_tree) {
	  list<Tree> result;
  	while (!first_list_of_arrows.is_this_the_end()) {
	  	auto enlarged_tree=tree; 
	  	enlarged_tree.add_node_with_arrows(first_list_of_arrows);
	  	result.splice(result.end(), add_nodes(enlarged_tree,nodes-1,first_list_of_arrows,std::forward<Closure>(is_valid_tree)));
		  ++first_list_of_arrows;
	  }
	  return result;
}

list<Tree> WaysToAddANode::add_intermediate_nodes(const Tree& tree, int nodes) {
  auto is_valid_tree =  [nodes] (const Tree& tree) {return tree.is_valid_enlargement_by(nodes);};
  return add_nodes(tree,nodes,first(),is_valid_tree);
}

list<Tree> WaysToAddANode::add_last_nodes(const Tree& tree, int nodes) {
  auto is_valid_tree =  [nodes] (const Tree& tree) {return tree.is_valid_completion_by(nodes);};
  return add_nodes(tree,nodes,first(),is_valid_tree);
}

ostream& operator<<(ostream& os, const WaysToAddANode& ways) {
		os<<"ways to add a node"<<endl;
		for (ListOfArrows list = ways.first();!list.is_this_the_end();++list) {
			os<<"- "; 
			for (int arrow: list) os<<arrow<<" ";
			os<<endl;
		}
		return os;
}	


template<typename T> bool is_permutation_of(vector<T> v, vector<T> w) {
	sort(v.begin(),v.end());
	sort(w.begin(),w.end());
	return equal(v.begin(),v.end(),w.begin(),w.end());
}


bool Tree::is_valid_enlargement_by(int nodes, bool verbose) const {
		if( verbose) {
			cout<<"checking validity of "<<endl; 
			cout<<*this;
			cout<<"nodes "<<nodes<<endl;
		}
		//check that all the original nodes have an incoming arrow.
		for (int i=0;i<no_nodes-nodes;++i)
			if (incoming_arrows[i]==0) {
				if (verbose) cout<<"not valid"<<endl;	
				return false;
			}
		if (verbose) cout<<"valid"<<endl;
		return true;
	}

bool Tree::is_valid_completion_by(int nodes, bool verbose) const {
		if( verbose) {
			cout<<"checking validity of "<<endl; 
			cout<<*this;
			cout<<"nodes "<<nodes<<endl;
		}
		//check that the original nodes all have an even, nonzero number of incoming arrows.
		for (int i=0;i<no_nodes-nodes;++i) {
			auto incoming_arrows_at_node=incoming_arrows[i] ;
			if (incoming_arrows_at_node==0 || incoming_arrows_at_node%2) {
				if (verbose) cout<<"not valid"<<endl;	
				return false;
			}
		}
		if (verbose) cout<<"valid"<<endl;
		return true;
	}







list<Tree> enlarge_tree(const Tree& tree, int nodes) {
//	for each new_node attach in all possible ways to an inner node
//	in such a way that each inner node is attached to something
	auto ways_to_add_a_node=tree.ways_to_add_a_node();
	return ways_to_add_a_node.add_intermediate_nodes(tree, nodes);
}

list<Tree> complete_tree(const Tree& tree, int nodes) {
	assert(nodes>0);
	auto ways_to_add_a_node=tree.ways_to_add_a_node();
	auto trees=ways_to_add_a_node.add_last_nodes(tree,nodes);
	return trees;
}



list<Tree> enlarge_incomplete_trees(int nodes_to_add, list<Tree>&& incomplete_trees) {
	list<Tree> result;
	for (auto& tree: incomplete_trees) 
		result.splice(result.end(),enlarge_tree(tree,nodes_to_add));
	return result;
}

list<Tree> incomplete_trees(const vector<int>& partition, int final_dimension) {
	assert(partition.size()>0);
	return (partition.size()==1) ?  
	  list<Tree>{Tree{partition[0],final_dimension}} :
	  enlarge_incomplete_trees(partition[0], incomplete_trees ({partition.begin()+1,partition.end()},final_dimension));
}

list<Tree> complete_incomplete_trees(int nodes_to_add, list<Tree>&& incomplete_trees) {
	list<Tree> result;
	nice_log<<incomplete_trees.size()<<"trees to complete"<<endl;
	remove_equivalent_trees(incomplete_trees);  //apparently useless (??)
	nice_log<<incomplete_trees.size()<<"trees to complete after removal of equivalents"<<endl;
  auto counter= nice_log.counter("completing tree",100);
	for (auto&& tree: move(incomplete_trees)) {
  	  nice_log<<counter;
  	  auto completions=complete_tree(tree,nodes_to_add);
  	  //nice_log<<completions.size()<<" completions"<<endl;
	  	result.splice(result.end(),completions);
  }
	return result;
}

list<Tree> trees(vector<int> partition) {
	assert(partition.size()>0);
	string label=get_label(partition);
	nice_log<<"computing trees with partition "<<label<<endl;
	int dimension=accumulate(partition.begin(),partition.end(),0);
	list<Tree> result = (partition.size()==1) ? 
	  list<Tree>{Tree{partition[0],dimension}} : 
	  complete_incomplete_trees(partition[0],incomplete_trees({partition.begin()+1,partition.end()},dimension));
	remove_equivalent_trees(result);
//	for (auto& tree : result)	tree.set_name(label);
	return result;
}




