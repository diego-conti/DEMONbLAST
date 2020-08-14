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
#ifndef FILTER_H
#define FILTER_H
#include "weightbasis.h"
#include "labeled_tree.h"
#include "includes.h"
#include "adinvariantobstruction.h"

//TODO derive from Options
class Filter {
	friend void DiagramDataOptions::adapt_to_filter(const Filter& filter);
  bool allow_nonnice=false;
  bool only_traceless_derivations_=false;
  int kernel_dimension_=-1;
  int cokernel_dimension_=-1;  
  bool only_nontrivial_automorphisms_=false;
  bool only_with_metric_=false;
	bool only_passing_obstruction_for_ad_invariant_metric_=false;
	tribool simple_nikolayevsky_=indeterminate;
public:
  void only_traceless_derivations()  {only_traceless_derivations_=true;}
	void only_nontrivial_automorphisms() {only_nontrivial_automorphisms_=true;}
	void only_with_metric() {only_with_metric_=true;}
	void only_passing_obstruction_for_ad_invariant_metric() {only_passing_obstruction_for_ad_invariant_metric_=true;}
  void N1N2N3() {allow_nonnice=true;}
  void simple_nikolayevsky(tribool choice) {simple_nikolayevsky_=choice;}
 	void only_kernel_MDelta_dimension(int kernel_dimension) {kernel_dimension_=kernel_dimension;}
 	void only_cokernel_MDelta_dimension(int cokernel_dimension) {cokernel_dimension_=cokernel_dimension;}
  bool trivial() const {
  	return !(only_traceless_derivations_||kernel_dimension_>=0||cokernel_dimension_>=0||only_nontrivial_automorphisms_ | only_with_metric_ | only_passing_obstruction_for_ad_invariant_metric_) && indeterminate(simple_nikolayevsky_);}
  bool meets(const LabeledTree& diagram, DiagramDataOptions options) const {
  	if (!allow_nonnice && !satisfies_formal_jacobi(diagram)) return false;
  	if (trivial()) return true;
    auto& properties = diagram.weight_basis(options).properties();
    if (only_traceless_derivations_ && !properties.are_all_derivations_traceless()) return false;
    if (kernel_dimension_>=0 && properties.dimension_kernel_M_Delta()!=kernel_dimension_) return false;
    if (cokernel_dimension_>=0 && properties.dimension_cokernel_M_Delta()!=cokernel_dimension_) return false;
    if (only_nontrivial_automorphisms_ && !properties.has_nontrivial_automorphisms()) return false;
    if (only_with_metric_ && !properties.potentially_admits_metrics()) return false;
    if (only_passing_obstruction_for_ad_invariant_metric_ && !passes_obstruction_for_ad_invariant_metric(diagram)) return false;
    if (simple_nikolayevsky_ && !properties.simple_nikolayevsky_derivation()) return false;
    if (!simple_nikolayevsky_ && properties.simple_nikolayevsky_derivation()) return false;
    return true;
  }
  
	bool require_only_nontrivial_automorphisms() const {return only_nontrivial_automorphisms_;}
	bool has_N1N2N3() const {return allow_nonnice;}

};

list<LabeledTree> nice_diagrams(vector<int> partition, const Filter& filter,DiagramDataOptions options);
void remove_trees(list<LabeledTree>& trees, const Filter& filter,DiagramDataOptions options);
#endif
