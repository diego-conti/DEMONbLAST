/*  Copyright (C) 2018-2023 by Diego Conti, diego.conti@unipi.it      
                                                                     
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

#include "includes.h"
#include "ricci.h"
#include "log.h"
#include "weightbasis.h"

using namespace Wedge;


class Derivations {
	ex generic_derivation_;
	set<ex,ex_is_less> remaining_equations_;
	VectorSpace<DifferentialForm> space_of_diagonal_derivations_,space_containing_offdiagonal_derivations_,space_contained_in_offdiagonal_derivations_;

	void compute_offdiag(const GL& gl, const LieGroup& G);
	void compute_diag(const GL& gl, const LieGroup& G);
public:
	Derivations(const GL& gl, const LieGroup& niceLieGroup);
	exvector remaining_equations() const {
		return {remaining_equations_.begin(), remaining_equations_.end()};
	}
	static exvector remaining_equations(const GL& gl, const LieGroup& niceLieGroup,ex generic_element); 
	bool always_a_derivation() const {
		return remaining_equations_.empty();
	}
	pair<int,int> dimension() const {
		return make_pair(space_of_diagonal_derivations_.Dimension()+space_contained_in_offdiagonal_derivations_.Dimension(),
			space_of_diagonal_derivations_.Dimension()+space_containing_offdiagonal_derivations_.Dimension());
	}
	VectorSpace<DifferentialForm> space_containing_offdiagonal_derivations() const {
		return space_containing_offdiagonal_derivations_;
	}
	VectorSpace<DifferentialForm> space_contained_in_offdiagonal_derivations() const {
		return space_contained_in_offdiagonal_derivations_;
	}
	VectorSpace<DifferentialForm> diagonal_derivations() const {
		return space_of_diagonal_derivations_;
	}
};

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
	Derivations derivations(const GL& gl) const;		//return the derivations as a subset of the Lie algebra gl
	string derivations() const;	
	
	exvector csquared(const WeightBasis& weight_basis) const;
	exvector c(const list<Weight>& weights) const;
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
