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
#include "filter.h"
#include "horizontal.h"

void remove_trees(list<LabeledTree>& trees, const Filter& filter,DiagramDataOptions options) {
  auto i=trees.begin();
  while (i!=trees.end()) 
      if (!filter.meets(*i,options)) 
      	 i=trees.erase(i); 
      else ++i;
}


list<LabeledTree> nice_diagrams(vector<int> partition, const Filter& filter,DiagramDataOptions options) {
  nice_log<<"nice diagrams of partition "<<horizontal(partition)<<endl;
  if (partition[0]<2) {
    nice_log<<"no trees"<<endl;
    return {};
  }
	list<LabeledTree> result;
	auto counter=nice_log.counter("Labeling tree ",50);
	for (auto& tree : trees(partition)) {
	  nice_log<<counter;
		LabeledTree labeled_tree{tree};
		auto labelings = labeled_tree.labelings();
		remove_trees(labelings, filter,options);
		remove_equivalent_trees_same_hash(labelings);
		result.splice(result.end(),labelings);
	}
	return result;
}

void DiagramDataOptions::adapt_to_filter(const Filter& filter) {
  	if(filter.require_only_nontrivial_automorphisms())
      set(DiagramDataOption::with_automorphisms);
    if (filter.only_with_metric() && !with_metrics()) throw runtime_error("Filter specifies metrics, but diagram data options do not specify the type of metric");
}

