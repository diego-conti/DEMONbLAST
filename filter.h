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
#include <vector>
using namespace std;

class Filter {
  bool allow_nonnice=false;
  bool no_coordinate_hyperplane_=false;
  bool only_traceless_derivations_=false;
  bool only_MDelta_surjective_=false;
  bool only_reachable_orthant_=false;
public:
  void no_coordinate_hyperplane()  {no_coordinate_hyperplane_=only_traceless_derivations_=true;}
  void only_traceless_derivations()  {only_traceless_derivations_=true;}
  void only_MDelta_surjective()  {only_MDelta_surjective_=true;}
  void only_reachable_orthant()  {only_reachable_orthant_=only_traceless_derivations_=true;}
  void N1N2N3() {allow_nonnice=true;}
  bool nontrivial() const {return no_coordinate_hyperplane_||only_traceless_derivations_||only_MDelta_surjective_||only_reachable_orthant_;}
  bool meets(const LabeledTree& diagram) const {
    return allow_nonnice || satisfies_formal_jacobi(diagram);
  }
  bool meets(const WeightBasis& weight_basis) const {
    auto& properties = weight_basis.properties();
    if (only_traceless_derivations_ && !properties.are_all_derivations_traceless()) return false;
    if (no_coordinate_hyperplane_ && properties.is_X_ijk_in_coordinate_hyperplane()) return false;
    if (only_MDelta_surjective_ && !properties.is_M_Delta_surjective()) return false;
    if (only_reachable_orthant_ && !properties.is_M_Delta_surjective()) throw invalid_argument("only-reachable-orthant not applicable when M_delta is not surjective");
    if (only_reachable_orthant_ && !static_cast<const DiagramPropertiesSurjectiveMDelta&>(properties).is_X_ijk_in_reachable_orthant()) return false;    
    return true;
  }
};

list<LabeledTree> nice_diagrams(vector<int> partition, const Filter& filter);

#endif
