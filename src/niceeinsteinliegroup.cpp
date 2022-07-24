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
#include "niceeinsteinliegroup.h"
#include "linearsolve.h"
#include "gauss.h"
#include "weightbasis.h"

list<NiceEinsteinLieGroup> NiceEinsteinLieGroup::from_coefficient_configuration (MetricCoefficientConfiguration&& configuration)
{
		list<NiceEinsteinLieGroup> result;
		insert_new_lie_group(result,configuration);
		while (configuration) {
			insert_new_lie_group(result,configuration);
			++configuration;
		}
		return result;

}



list<NiceEinsteinLieGroup> NiceEinsteinLieGroup::from_weight_basis(const WeightBasisAndProperties& weight_basis,MetricType metric_type) {
		MetricCoefficientConfiguration configuration{weight_basis,metric_type};
    list<NiceEinsteinLieGroup> result=from_coefficient_configuration(move(configuration));
    return result;
}

NiceEinsteinLieGroup::NiceEinsteinLieGroup(const MetricCoefficientConfiguration& configuration) : LieGroupsFromDiagram(configuration.lie_algebra_dimension()) {
		ExVector de(configuration.lie_algebra_dimension());
		int i=0;
		for (WeightAndValue weight : configuration.weights()) 
			de[weight.node_out]+=weight.value*e()[weight.node_in1]*e()[weight.node_in2];
		for (int i=1;i<=de.size();++i)
			Declare_d(e(i), de(i));
}



void NiceEinsteinLieGroup::insert_new_lie_group(list<NiceEinsteinLieGroup>& out_list, const MetricCoefficientConfiguration& configuration) {
	NiceEinsteinLieGroup nice_lie_group{configuration};
	if (nice_lie_group.solve_linear_ddzero()) {
	  nice_log<<horizontal(nice_lie_group.StructureConstants())<<endl;
	  for (auto x: nice_lie_group.e()) {
	    auto ddx=nice_lie_group.d(nice_lie_group.d(x));
	    if (!ddx.is_zero()) nice_log<<"d^2!=0 in general, "<<ddx<<endl;
	  }
		out_list.push_back(move(nice_lie_group));
	}
}


