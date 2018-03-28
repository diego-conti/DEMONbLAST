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
#ifndef LABELED_TREE_H
#define LABELED_TREE_H
#include "tree.h"
#include <memory>

using namespace std;
class WeightBasis;

class LabeledTree : public TreeBase<LabeledArrow> {
public:
	using TreeBase<LabeledArrow>::TreeBase;
	LabeledTree(const Tree& tree);
	LabeledTree(const LabeledTree& labeled_tree);
	list<LabeledTree> labelings() const;
	list<Weight> weights() const {
		list<Weight> result;
		for (auto arrow : arrows()) 
			if (arrow.node_in<arrow.label) 
				result.push_back({arrow.node_in,arrow.label,arrow.node_out});
		return result;
	}
	const WeightBasis& weight_basis() const;
	virtual ~LabeledTree(); 
private:
	friend struct TestLabeledTree;
	friend class WeightedTree;

  mutable unique_ptr<WeightBasis> weight_basis_;
	LabeledTree with_added_label(int to_node, int from_1, int from_2) const;
	list<LabeledTree> labelings_with_bracket(int to_node, int from_1, int from_2) const;
  list<LabeledTree> labelings_with_hint(int hint_node) const;
	vector<int> unlabeled_arrows_pointing_to(int node) const;
  vector<int> number_of_unlabeled_incoming_arrows() const;
	int node_with_less_unlabeled_incoming_arrows() const;
	LabeledTree() = default;
	void add_arrow(LabeledArrow&& arrow) = delete;
};

bool satisfies_formal_jacobi(const LabeledTree& tree);
bool has_double_incoming_arrow_with_multiplicity_three(const LabeledTree& tree);

#endif
