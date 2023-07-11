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
#include "liegroupsfromdiagram.h"
#include "linearsolve.h"
#include "gauss.h"
#include "weightbasis.h"
#include "horizontal.h"
#include "linearinequalities.h"

using namespace Wedge;
using namespace std;

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
    for (ex eq : eqns) eqns2.append(eq.expand().numer());
		return impose_polynomial_eqns<StructureConstant>(*this, move(eqns2), is_solution_in_positive_orthant);
}


exvector LieGroupsFromDiagram::csquared(const WeightBasis& weight_basis) const {
		exvector result;
		for (auto& weight: weight_basis.weights_and_coefficients()) 
			result.push_back(pow(c_ijk(weight.node_in1,weight.node_in2,weight.node_out),2));
		return result;
}
exvector LieGroupsFromDiagram::c(const list<Weight>& weights) const {
		exvector result;
		for (auto& weight: weights) 
			result.push_back(-c_ijk(weight.node_in1,weight.node_in2,weight.node_out));
		return result;
}



class ContainerOfEquations {
	exvector eqns;
public:
	template<typename Container>
	ContainerOfEquations(Container&& container) {
		for (auto eq: container) insert(eq);
	}
	void insert(ex eq) {
		eq=eq.expand().numer();
		if (eq.is_zero()) return;
		if (find(eqns.begin(),eqns.end(),eq)!=eqns.end()) return;
		if (find(eqns.begin(),eqns.end(),-eq)!=eqns.end()) return;
		eqns.push_back(eq);
	}
	auto begin() const {return eqns.begin();}
	auto end() const {return eqns.end();}
	bool empty() const {return eqns.empty();}
	string to_string() const {
	  stringstream ss;
	  ss<<latex;
	  for (auto eq:eqns) ss<<eq<<"=0; ";
 	  return ss.str();	
	}
};

string to_string(const LieGroupHasParameters<true>& G) {
	stringstream s;
	s<<latex;
	G.canonical_print(s);
	set<ex,ex_is_less> eqns;
	G.GetEquations_ddZero(eqns);
	auto eqns2=ContainerOfEquations{eqns};
	if (!eqns2.empty()) s<<" d^2=0 when "<<eqns2.to_string();
	return s.str();
}

bool LieGroupsFromDiagram::is_dd_nonzero() const {
		lst eqns;
		GetEquations_ddZero(eqns);
		for (auto eq: eqns) if (is_a<numeric>(eq) && !eq.is_zero()) return true;
		return false;
}

ex Xbracket(const LieGroup& G, const GLRepresentation<VectorField>& V, ex A, ex X, ex Y) {
	ex Ax=V.Action<VectorField>(A,X);
	ex Ay=V.Action<VectorField>(A,Y);
	ex Axy=V.Action<VectorField>(A,G.LieBracket(X,Y));
	return G.LieBracket(Ax,Y)+G.LieBracket(X,Ay)-Axy;
}

exvector Xbrackets(const LieGroup& G, const GLRepresentation<VectorField>& V, ex A) {
		exvector Xbrackets;
		for (int i=1;i<=G.Dimension();++i)
		for (int j=i+1;j<=G.Dimension();++j) 
			Xbrackets.push_back(Xbracket(G,V,A,G.e(i),G.e(j)));				
		return Xbrackets;
}

VectorSpace<DifferentialForm> offdiagonal_elements(const GL& gl) {
	exvector basis;
	for (int i=1;i<=gl.n();++i)
	for (int j=i+1;j<=gl.n();++j) {
		basis.push_back(gl.A(i,j));
		basis.push_back(gl.A(j,i));
	}
	return {basis.begin(),basis.end()};
}

VectorSpace<DifferentialForm> diagonal_elements(const GL& gl) {
	exvector basis;
	for (int i=1;i<=gl.n();++i)
		basis.push_back(gl.A(i,i));
	return {basis.begin(),basis.end()};
}


auto SubspaceFromSolutions(const VectorSpace<DifferentialForm>& V, const lst& sol) {
	list<ex> equations;
	for (auto eq: sol) equations.push_back(eq.lhs()-eq.rhs());
	return V.SubspaceFromEquations(equations.begin(),equations.end());
}


exvector Derivations::remaining_equations(const GL& gl, const LieGroup& G,ex generic_element) {
	auto X=Xbrackets(G,GLRepresentation<VectorField>(&gl,G.e()),generic_element);
	exvector eqns;
	GetCoefficients<VectorField>(eqns,X);
	return eqns;
}

void Derivations::compute_offdiag(const GL& gl, const LieGroup& G) {
	auto V=offdiagonal_elements(gl);
	auto X=Xbrackets(G,GLRepresentation<VectorField>(&gl,G.e()),V.GenericElement());
	lst eqns;
	GetCoefficients<VectorField>(eqns,X);
	Wedge::linear_impl::LinearEquationsWithParameters<VectorSpace<DifferentialForm>::Coordinate,StructureConstant> equations(eqns);
	while (equations.eliminate_linear_equations()) ;
	auto sol=equations.solution();
	space_containing_offdiagonal_derivations_=SubspaceFromSolutions(V,sol);
	sol=equations.always_solution();
	space_contained_in_offdiagonal_derivations_=SubspaceFromSolutions(V,sol);
	X=Xbrackets(G,GLRepresentation<VectorField>(&gl,G.e()),space_containing_offdiagonal_derivations_.GenericElement());
	GetCoefficients<VectorField>(remaining_equations_,X);
		remaining_equations_.erase(0);					
}

void Derivations::compute_diag(const GL& gl, const LieGroup& G) {
	auto V=diagonal_elements(gl);
	auto X=Xbrackets(G,GLRepresentation<VectorField>(&gl,G.e()),V.GenericElement());
	lst eqns;
	GetCoefficients<VectorField>(eqns,X);
	space_of_diagonal_derivations_=V.SubspaceFromEquations(eqns.begin(),eqns.end());
}

Derivations::Derivations(const GL& gl, const LieGroup& G)  {
	compute_diag(gl,G);
	compute_offdiag(gl,G);
}

Derivations LieGroupsFromDiagram::derivations(const GL& gl) const {
		return Derivations{gl,*this};
}

string LieGroupsFromDiagram::derivations() const {
	GL Gl(Dimension());
	auto Der=derivations(Gl);
		stringstream s;
		pair<int,int> dim=Der.dimension();
		if (Der.always_a_derivation())
			s<<"dim Der(g)="<<dim.first<<endl;		
		else 
			s<<dim.first<<"<= dim Der(g)<="<<dim.second<<", derivation only if "<<horizontal(Der.remaining_equations())<<endl;
		matrix generic_offdiag_derivation=Gl.glToMatrix(Der.space_containing_offdiagonal_derivations().GenericElement());
		if (generic_offdiag_derivation.pow(Dimension()).is_zero_matrix()) s<<" offdiag derivations are nilpotent"<<endl;
		else s<<latex<<" offdiag derivation are not nilpotent: "<<generic_offdiag_derivation<<endl;
		return s.str();
};

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


