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
#include "tree.h"
#include "labeled_tree.h"

class DoubleArrowsTree : public TreeBase<LabeledDoubleArrow> {
public:
  DoubleArrowsTree(const LabeledTree& tree) : TreeBase(tree.number_of_nodes(),tree.number_of_nodes()) {
    auto weights=tree.weights();
  	for (auto i=weights.begin();i!=weights.end();++i) {
			auto j=i;
			while (++j!=weights.end()) {
				add_concatenated_arrow(*i,*j);	
				add_concatenated_arrow(*j,*i);	
			}
		}  
  };
private:
  void add_arrow_if_no_repetition(int from, int to, pair<int,int> label) {
    if (from!=label.first && from!=label.second) add_arrow({from,to,label});
  }
  void add_concatenated_arrow(Weight i, Weight j) {
    if (i.node_out == j.node_in1) add_arrow_if_no_repetition(j.node_in2, j.node_out,make_pair(i.node_in1,i.node_in2));
    else if (i.node_out == j.node_in2) add_arrow_if_no_repetition(j.node_in1, j.node_out,make_pair(i.node_in1,i.node_in2));
  }
};



