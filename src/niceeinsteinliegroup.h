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
#ifndef NICE_EINSTEIN_LIE_GROUP_H
#define NICE_EINSTEIN_LIE_GROUP_H
#include "tree.h"
#include "labeled_tree.h"
#include "xginac.h"
#include "weightbasis.h"
#include "liegroupsfromdiagram.h"
#include "horizontal.h"
#include "coefficientconfiguration.h"

using namespace Wedge;


class NiceEinsteinLieGroup : public LieGroupsFromDiagram {
	NiceEinsteinLieGroup(const MetricCoefficientConfiguration& configuration);	
  static	list<NiceEinsteinLieGroup> from_coefficient_configuration (MetricCoefficientConfiguration&& configuration);
 	static void insert_new_lie_group(list<NiceEinsteinLieGroup>& out_list, const MetricCoefficientConfiguration& configuration);
public:
	static list<NiceEinsteinLieGroup> from_weight_basis(const WeightBasisAndProperties& weight_basis, MetricType);
};




#endif
