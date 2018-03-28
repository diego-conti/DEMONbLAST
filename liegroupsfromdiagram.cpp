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
#include "liegroupsfromdiagram.h"
#include "linearsolve.h"
#include "gauss.h"
#include "weightbasis.h"
#include "wedge/liesubgroup.h"
#include "horizontal.h"
#include "linearinequalities.h"


using namespace Wedge;

bool LieGroupsFromDiagram::solve_linear_ddzero() {
		static auto is_solution_in_positive_orthant = [] (ex sol) {
			exvector rhs;
			transform(sol.begin(),sol.end(),back_inserter(rhs),[] (ex solution) {return solution.rhs();});
			bool positive=LinearInequalities {rhs.begin(),rhs.end(), StructureConstant{}}.has_solution();
			return positive;
		};
		lst eqns;
		GetEquations_ddZero(eqns);
    lst eqns2;
    for (ex eq : eqns) eqns2.append(eq.expand());
		return impose_polynomial_eqns<StructureConstant>(*this, move(eqns), is_solution_in_positive_orthant);
}



template<typename IteratorBegin, typename IteratorEnd>
string ddzero(IteratorBegin begin, IteratorEnd end) {
  stringstream ss;
  while (begin!=end) ss<<*begin++<<"=0; ";
  return ss.str();
}

string to_string(const LieGroupHasParameters<true>& G) {
  string structure_constants =horizontal(G.StructureConstants());
  set<ex,ex_is_less> eqns;
	G.GetEquations_ddZero(eqns);
	eqns.erase(0);
	return eqns.empty()?  structure_constants :  structure_constants+" d^2=0 when "+ddzero(eqns.begin(),eqns.end());
}



bool LieGroupsFromDiagram::is_dd_nonzero() const {
		lst eqns;
		GetEquations_ddZero(eqns);
		for (auto eq: eqns) if (is_a<numeric>(eq) && !eq.is_zero()) return true;
		return false;
}

AbstractLieSubgroup<true> inverted_structure_constants(const LieGroupsFromDiagram& G) {
  ExVector frame;
  for (auto X : G.e()) frame.insert(frame.begin(),-X);
  return {G,frame};
}
AbstractLieSubgroup<true> change_basis(const LieGroupsFromDiagram& G, const vector<int>& sigma) {
  ExVector frame;
  for (int sigma_i : sigma) frame.push_back(G.e(sigma_i+1));
  return {G,frame};
}


exvector generic_metric(int dimension) {
  exvector result;
  for (int i=1;i<=dimension;++i)
    result.push_back(realsymbol("epsilon"+to_string(i)));
  return result;
}
ex identity_matrix(int order) {
  matrix I(order,order);
   for (int i=0;i<order;++i) I(i,i)=1;
  return I;
}

template<typename Container, typename Separator, typename... Separators> 
string to_string(Container&& container, Separator first_separator, Separators... other_separators) {
  if (container.empty()) return {};
  string result;
  auto i=container.begin();
  result=to_string(*i,other_separators...);
  while (++i!=container.end()) result+=first_separator+to_string(*i,other_separators...);
  return result;
}

    
string nontrivial_automorphisms_to_string(list<vector<int>>&& automorphisms) {
  for (auto& automorphism : automorphisms) 
    for (int& element : automorphism) ++element;
  if (automorphisms.empty()) return {};
  else return "Nontrivial automorphisms:\n" +   to_string(automorphisms,"\n",",") + "\n";
}


