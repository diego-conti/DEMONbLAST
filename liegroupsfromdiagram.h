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
#ifndef LIE_GROUPS_FROM_DIAGRAM_H
#define LIE_GROUPS_FROM_DIAGRAM_H
#include <wedge/liegroup.h>
#include <wedge/liesubgroup.h>
#include <wedge/gl.h>
#include "ricci.h"
#include "log.h"
#include "weightbasis.h"

using namespace Wedge;


class LieGroupsFromDiagram : public LieGroupHasParameters<true>, public ConcreteManifold, public virtual Has_dTable {
protected:
	bool solve_linear_ddzero();
	ex c_ijk(int i, int j, int k) const {
		return Hook(e()[j]*e()[i],d(e()[k]));
	}
public:
	using ConcreteManifold::ConcreteManifold;
	using Has_dTable::Declare_d;
	bool is_dd_nonzero() const;
	void DeclareConditions(const lst& list_of_equations) override {
		for (exmap::const_iterator i=dTable().begin();i!=dTable().end();i++)					
			Has_dTable::Declare_d(i->first,i->second.subs(list_of_equations));
	}
	VectorSpace<DifferentialForm> derivations(const GL& gl) const;		//return the derivations as a subset of the Lie algebra gl
	string derivations() const;	
	
	exvector csquared(const WeightBasis& weight_basis) const;
};


string to_string(const LieGroupHasParameters<true>& G);

AbstractLieSubgroup<true> change_basis(const LieGroupsFromDiagram& G, const vector<int>& sigma);
AbstractLieSubgroup<true> inverted_structure_constants(const LieGroupsFromDiagram& G);


exvector generic_metric(int dimension); 
ex identity_matrix(int order);

template<typename NiceGroupClass>
void test_einstein_condition(const list<NiceGroupClass>& groups, const exvector& diagonal_metric) {
  for (auto group: groups) {    
    int dimension=group.Dimension();    
    assert(dimension==diagonal_metric.size());
    
    
    auto ricci=ricci_operator(group,diagonal_metric);
    ex s=ricci.op(0);
    if (!ex_to<matrix>((ricci-s*identity_matrix(dimension)).evalm()).is_zero_matrix()) {
      nice_log<<"ERROR: not Einstein"<<endl<<ricci<<endl;
      nice_log<<ricci_operator(group,generic_metric(dimension));
    }
    else {
      nice_log<<"OK : ricci = "<<ricci<<endl;
    }
  }
}


template<typename NiceGroupClass>
void test_einstein_condition(const list<NiceGroupClass>& groups) {
 for (auto group: groups) {
    nice_log<<to_string(group)<<endl;
    int dimension=group.Dimension();
    nice_log<<"ricci"<<endl<<ricci_operator(group,generic_metric(dimension))<<endl;
   }
}

string nontrivial_automorphisms_to_string(list<vector<int>>&& automorphisms);


#endif
