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
#ifndef NICE_H
#define NICE_H

#include "arrow.h"
#include "permutations.h"
#include "log.h"
#include "horizontal.h"


class Tree;

namespace Position {
  class begin_t {};
  class end_t {};
  static begin_t begin {};
  static end_t end {};
}

class ListOfArrows {
  using BooleanType = char; 
  struct ListOfArrowsConstIterator {
    const ListOfArrows& list_of_arrows;
    int index;
    ListOfArrowsConstIterator(const ListOfArrows& list, Position::begin_t) : list_of_arrows{list}, index{0} {
      while (!list_of_arrows.at(index) && index<list.size()) ++index;
    }
    ListOfArrowsConstIterator(const ListOfArrows& list,  Position::end_t) : list_of_arrows{list}, index{list.size()} {}
    int operator*() const {return index;}  
    ListOfArrowsConstIterator& operator++();
    bool operator!=(ListOfArrowsConstIterator iter) {return index!=iter.index;}
    void print(ostream& os) {cout<<"iterator points to "<<index<<endl;}
  };
public:
  ListOfArrows(int nodes, Position::begin_t) : arrow_present(nodes+1) {}
  ListOfArrows(int nodes, Position::end_t) : arrow_present(nodes+1) {arrow_present.back()=true;}
  ListOfArrows& operator++();
  ListOfArrowsConstIterator begin() const {return ListOfArrowsConstIterator{*this,Position::begin};}
  ListOfArrowsConstIterator end() const {return ListOfArrowsConstIterator{*this,Position::end};}  
  bool is_this_the_end() const {return arrow_present.back();}
  bool at(int index) const {return arrow_present[index];}
  int size() const {return arrow_present.size()-1;}
  bool operator!=(const ListOfArrows& list_of_arrows) {return arrow_present!=list_of_arrows.arrow_present;}
private:
//contains one extra element which is set to true to represent the list_of_arrows "past the end" when iterating through lists of arrows
  vector<BooleanType> arrow_present; 
};



class WaysToAddANode {
public:
  WaysToAddANode(int number_of_nodes) : nodes{number_of_nodes} {}
	list<Tree> add_intermediate_nodes(const Tree& tree, int nodes);
	list<Tree> add_last_nodes(const Tree& tree, int nodes);
	ListOfArrows first() const {return ListOfArrows{nodes,Position::begin};}
private:
  int nodes;
};

ostream& operator<<(ostream& os, const WaysToAddANode& ways);

class FixedLengthVectorInt : private vector<int> {
public:
  using vector<int>::vector;
  using vector<int>::operator[];
  vector<int> to_vector(int size) const {return {begin(),begin()+size};}
  bool has_at_least_capacity(int size) const {
    return vector<int>::size()>=size;
  }
};

class SetOfNodes {
public:
	SetOfNodes(int nodes, int final_size) : no_nodes{nodes},  nodes_hash(final_size) {}
	SetOfNodes()=default;
	string name() const {
	  if (name_.empty()) return std::to_string(number_)+ "#"+std::to_string(tree_hash);
	  else  return std::to_string(number_)+ "#"+std::to_string(tree_hash)+ " : "+name_;
	}
	void set_name(string new_name) {name_=move(new_name);}
	void add_number_to_name(int n) {number_=n;}
	int number_of_nodes() const {return no_nodes;}
	bool is_hash_computed() const {return tree_hash!=HASH_NOT_COMPUTED;}
protected:
	static const int HASH_NOT_COMPUTED=0;
	int no_nodes=0;
	mutable int tree_hash=HASH_NOT_COMPUTED;
	mutable FixedLengthVectorInt nodes_hash;
	bool capacity_exceeded() {return !nodes_hash.has_at_least_capacity(no_nodes);}
private:
  int number_=0;
	string name_;
};


template<typename ArrowType> class TreeBase : public SetOfNodes {
private:
	list<ArrowType> arrows_in_tree;
	void compute_hash_and_cache_result() const;
protected:
	static const int NO_NODE=-1;
	void invalidate() {tree_hash=HASH_NOT_COMPUTED;}
	bool matches(const TreeBase& tree, const vector<int>& permutation) const;
  void add_arrow(const ArrowType& arrow) {
    arrows_in_tree.push_back(arrow);
		invalidate();
  }
  template<typename Arrow>
  TreeBase(const SetOfNodes& tree, const list<Arrow>& arrows) : SetOfNodes(tree),  arrows_in_tree{arrows.begin(),arrows.end()} {
    assert(tree_hash!=HASH_NOT_COMPUTED);
  }
  TreeBase(const SetOfNodes& tree, const list<ArrowType>& arrows) : SetOfNodes(tree),  arrows_in_tree{arrows} {
    assert(tree_hash!=HASH_NOT_COMPUTED);
  }
  TreeBase(const SetOfNodes& tree, list<ArrowType>&& arrows) : SetOfNodes(tree),  arrows_in_tree{move(arrows)} {
    assert(tree_hash!=HASH_NOT_COMPUTED);
  }
public:
	TreeBase(const TreeBase&)=default;
	TreeBase(TreeBase&&)=default;
	TreeBase() = default;
	using SetOfNodes::SetOfNodes;
  const list<ArrowType>& arrows() const {return arrows_in_tree;} 
	bool is_equivalent_to(const TreeBase& tree) const;
	list<vector<int>> nontrivial_automorphisms() const;
	vector<int> node_hash() const {
	  if (tree_hash==HASH_NOT_COMPUTED) compute_hash_and_cache_result();
	  return nodes_hash.to_vector(number_of_nodes());
	}
	int hash() const {
	  if (tree_hash==HASH_NOT_COMPUTED) compute_hash_and_cache_result();
    return tree_hash;	  
	}
	string to_dot_string(string extra_data={}) const;
	void invert_nodes() {
    for (int i=0;2*i<number_of_nodes();++i)
      swap(nodes_hash[i],nodes_hash[number_of_nodes()-1-i]);
	  for (auto& arrow: arrows_in_tree)
	    arrow.invert(number_of_nodes());	 
   }
  template<typename Compare> 
	void sort_arrows(Compare&& compare) {
		arrows_in_tree.sort(compare);
	}
};

template<typename ArrowType>
ostream& operator<< (ostream& os, const TreeBase<ArrowType>& tree) {
  return os<<tree.to_dot_string();
}

class Tree : public TreeBase<Arrow> {
	friend class LabeledTree;
	FixedLengthVectorInt incoming_arrows;
public:
  Tree()=default;
  Tree(const Tree& tree) :  TreeBase<Arrow>(tree), incoming_arrows{tree.incoming_arrows} {}
  Tree(Tree&&)=default;
	Tree(int nodes, int final_size) : TreeBase<Arrow>(nodes,final_size), incoming_arrows(final_size) {}
	bool is_valid_enlargement_by(int nodes, bool verbose=false) const;
	bool is_valid_completion_by(int nodes, bool verbose=false) const;
	WaysToAddANode ways_to_add_a_node() const {return {no_nodes};}

	void add_node_with_arrows(const ListOfArrows& destinations) {
		int last_node=no_nodes++;
		assert(!capacity_exceeded());		
		for (auto destination: destinations) {
			add_arrow(Arrow{last_node,destination});
		  ++incoming_arrows[destination];
		}
		invalidate();
	}
};

template<typename Tree> 
void remove_equivalent_trees(list<Tree>& list_of_trees);

template<typename Tree> 
void remove_equivalent_trees_same_hash(list<Tree>& list_of_trees);

struct Weight {
	int node_in1, node_in2, node_out;
};

inline ostream& operator<<(ostream& os, Weight weight) {
	return os<<"{"<<weight.node_in1+1<<","<<weight.node_in2+1<<"}->"<<weight.node_out+1;
}

list<Tree> trees(vector<int> partition);

#include "tree.hpp"
#endif
