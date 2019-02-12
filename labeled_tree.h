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
#include "options.h"
#include <memory>

using namespace std;
class WeightBasisAndProperties;


enum class DiagramDataOption : unsigned int {
  dflt=0,
  no_data=0,
  with_diagonal_ricci_flat_metrics=1,
  with_sigma_compatible_ricci_flat_metrics=2,
  with_diagonal_nilsoliton_metrics=4,
  with_im_delta2=8,
  analyze_diagram=16,
  with_matrix_data=32,
  with_antidiagonal_ricci_flat_sigma=64
};

struct DiagramDataOptions : Options<DiagramDataOption> {
  bool with_diagonal_ricci_flat_metrics() const {return has(DiagramDataOption::with_diagonal_ricci_flat_metrics);}
	bool with_diagonal_nilsoliton_metrics() const {return has(DiagramDataOption::with_diagonal_nilsoliton_metrics);}
  bool with_sigma_compatible_ricci_flat_metrics() const {return has(DiagramDataOption::with_sigma_compatible_ricci_flat_metrics);}
  bool with_im_delta2() const {return has(DiagramDataOption::with_im_delta2);}
  bool analyze_diagram() const {return has(DiagramDataOption::analyze_diagram);}
	bool with_data() const {return *this!=DiagramDataOption::no_data;}
	bool with_matrix_data() const {return has(DiagramDataOption::with_matrix_data);}
	bool with_antidiagonal_ricci_flat_sigma() const {return has(DiagramDataOption::with_antidiagonal_ricci_flat_sigma);}
	bool with_metrics() const {return with_diagonal_ricci_flat_metrics() || with_diagonal_nilsoliton_metrics() || with_sigma_compatible_ricci_flat_metrics();}
};


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
	const WeightBasisAndProperties& weight_basis(DiagramDataOptions diagram_data_options) const;
	void canonicalize_order_increasing();
	void canonicalize_order_decreasing();
	virtual ~LabeledTree(); 	
	static LabeledTree from_stream(istream& s);
	string as_string() const;
private:
	friend struct TestLabeledTree;
	friend class WeightedTree;

  mutable unique_ptr<WeightBasisAndProperties> weight_basis_;
	LabeledTree with_added_label(int to_node, int from_1, int from_2) const;
	list<LabeledTree> labelings_with_bracket(int to_node, int from_1, int from_2) const;
  list<LabeledTree> labelings_with_hint(int hint_node) const;
	vector<int> unlabeled_arrows_pointing_to(int node) const;
  vector<int> number_of_unlabeled_incoming_arrows() const;
	int node_with_less_unlabeled_incoming_arrows() const;
	LabeledTree() = default;
	using TreeBase<LabeledArrow>::add_arrow;
};

bool satisfies_formal_jacobi(const LabeledTree& tree);
bool has_double_incoming_arrow_with_multiplicity_three(const LabeledTree& tree);

#endif
