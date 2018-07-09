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
#include "niceliegroup.h"
#include "linearsolve.h"
#include "gauss.h"
#include "weightbasis.h"


list<NiceLieGroup> NiceLieGroup::from_coefficient_configuration(CoefficientConfiguration&& configuration) {
		list<NiceLieGroup> result;
		while (configuration) {
			insert_new_lie_group(result,configuration);
			++configuration;
		}	
		return result;
}

list<NiceLieGroup> NiceLieGroup::from_weight_basis(const WeightBasis& weight_basis) {
	CoefficientConfigurationWithoutRedundantParameter configuration{weight_basis};
	return from_coefficient_configuration(std::move(configuration));
}

NiceLieGroup::NiceLieGroup(const CoefficientConfiguration& configuration) : LieGroupsFromDiagram(configuration.lie_algebra_dimension()) {
		ExVector de(configuration.lie_algebra_dimension());
		int i=0;
		for (WeightAndValue weight : configuration.weights()) {
			de[weight.node_out]+=weight.value*e()[weight.node_in1]*e()[weight.node_in2];
			//if (einstein.X_ijk!=matrix{}) X_ijk.push_back(einstein.X_ijk(i++,0)/(weight.value*weight.value));
		}
		for (int i=1;i<=de.size();++i)
			Declare_d(e(i), de(i));
}
	
exvector de(const LieGroup& G) {
	exvector de;
	for (ex e : G.e()) de.push_back(G.d(e));
	return de;
}

lst substitutions(const exvector& variables, const exvector& values) {
	lst subs;
	assert(variables.size()==values.size());
	for (int i=0;i<variables.size();++i)
		subs.append(variables[i]==values[i]);
	return subs;
}

exvector subs_in_vector(exvector vector, const exvector& variables, const exvector& values) {
	lst subs=substitutions(variables,values);
	for (auto& x: vector)
		x=x.subs(subs).expand();
	return vector;
}

bool matches(const lst& substitutions, const exvector& lhs, const exvector& rhs) {
	assert(lhs.size()==rhs.size());
	for (int i=0;i<lhs.size();++i)
		if (!(lhs[i].subs(substitutions)-rhs[i]).expand().is_zero()) {
			nice_log<<"mismatch "<<horizontal(lhs)<<" != "<<horizontal(rhs)<<endl;
			nice_log<<"mismatch "<<(lhs[i].subs(substitutions)-rhs[i]).expand()<<endl;
			return false;
		}
	return true;
}

bool equivalent(const NiceLieGroup& G, const NiceLieGroup& H) {
	exvector constants;
	exvector de_G=de(G), de_H=de(H);	
	lst two_forms_G_to_H=substitutions(G.pForms(2).e(),H.pForms(2).e());
	if (matches(two_forms_G_to_H,de_G,de_H)) {
			nice_log<<"eliminating "<<to_string(H)<<", identical to "<<to_string(G)<<endl;
			return true;
	}
	GetSymbols<StructureConstant>(constants, de_G.begin(),de_G.end());
	nice_log<<"constants "<<constants<<endl;
	SignConfiguration signs_of_parameters{constants.size()};
	while (signs_of_parameters.has_next()) {
		++signs_of_parameters;
		nice_log<<"testing "<<signs_of_parameters<<endl;
		auto de_G_eq=subs_in_vector(de_G,constants,signs_of_parameters.multiply(constants));
		if (matches(two_forms_G_to_H,de_G_eq,de_H)) {
			nice_log<<"eliminating "<<to_string(H)<<", equivalent to "<<to_string(G)<<endl;
			return true;			
		}
	} 
	
	nice_log<<"mismatch "<<to_string(H)<<", not equivalent to "<<to_string(G)<<endl;
	return false;
}	


void NiceLieGroup::insert_new_lie_group(list<NiceLieGroup>& out_list, const CoefficientConfiguration& configuration) {
	NiceLieGroup nice_lie_group{configuration};
  nice_log<<to_string(nice_lie_group)<<endl;
	if (nice_lie_group.solve_linear_ddzero()) {
    nice_log<<"solve_linear_ddzero passed"<<endl;
		if (all_of(out_list.begin(),out_list.end(),[&nice_lie_group] (const NiceLieGroup& already_in_list) {
			return !equivalent(nice_lie_group,already_in_list);
		}))	out_list.push_back(move(nice_lie_group));
	}		
	else nice_log<<"solve_linear_ddzero failed"<<endl;
}


