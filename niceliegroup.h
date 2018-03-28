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
#ifndef NICE_LIE_GROUP_H
#define NICE_LIE_GROUP_H
#include "tree.h"
#include "labeled_tree.h"
#include <wedge/liegroup.h>
#include "liegroupsfromdiagram.h"
#include "xginac.h"
#include "weightbasis.h"

using namespace std;
using namespace Wedge;

class CoefficientConfiguration;

class NiceLieGroup : public LieGroupsFromDiagram {
	NiceLieGroup(const CoefficientConfigurationWithoutRedundantParameter& configuration);	
	static void insert_new_lie_group(list<NiceLieGroup>& out_list, const CoefficientConfigurationWithoutRedundantParameter& configuration);
  void DeclareConditions(const lst& list_of_equations) override {
    LieGroupsFromDiagram::DeclareConditions(list_of_equations);
	  for (auto& X : X_ijk) X=X.subs(list_of_equations);
  }
public:
	static list<NiceLieGroup> from_weight_basis(const WeightBasis& weight_basis);
  exvector X_ijk;
};



#endif
