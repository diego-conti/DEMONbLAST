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
#ifndef LINEARSOLVE_H
#define LINEARSOLVE_H

#include "includes.h"
#include "log.h"
namespace Wedge {
namespace linear_impl {

template<typename IsAdmissibleSolution>
lst solve_linear(const lst& eqns, const lst& unknowns,const IsAdmissibleSolution& isAdmissibleSolution) {
	ex sol=DefaultLinAlgAlgorithms::lsolve(eqns,unknowns);
	if (sol==lst() || !isAdmissibleSolution(sol)) return {};
	assert(is_a<lst>(sol));
	return ex_to<lst>(sol);
}

inline bool is_trivial_solution(const lst& sol) {
	return all_of(sol.begin(),sol.end(),[](ex eq) {return eq.rhs()==eq.lhs();});
}

template<class Variable> lst get_variables(const lst& eqns) {
	list<ex> variables;
	GetSymbols<Variable>(variables,eqns.begin(),eqns.end());
  return variables;
}


//represents polynomial equations possibly involving parameters
template<class Variable>
class AbstractPolynomialEquations {
  void update_solution(lst solution) {
    nice_log<<"solution: "<<solution<<endl;
    if (solution==lst{}) sol.remove_all();
    for (int i=0;i<sol.nops();++i)
			sol.let_op(i)=sol.op(i).lhs()==sol.op(i).rhs().subs(solution);
    for (int i=0;i<equations.nops();++i)
			equations.let_op(i)=equations.op(i).subs(sol).expand();
  }
  static lst expand(const lst& eqns) {
    lst result;
    for (auto eq: eqns) result.append(eq.expand());
    return result;
  }
  static lst expand(lst&& eqns) {
  	ex subs = abs(-wild()) == abs(wild());
    for (int i=0;i<eqns.nops();++i)
      eqns.let_op(i)=eqns.op(i).expand().subs(subs).expand();
    return std::move(eqns);
  }
protected:
  lst equations;
  lst variables;
  lst sol;
  virtual lst linear_equations() const =0;
<<<<<<< HEAD:src/linearsolve.h
  template<typename T>
  AbstractPolynomialEquations(std::initializer_list<T>)=delete;	//disable construction with {} to avoid ambiguities
public:
  template<typename ListOfEquations>
  AbstractPolynomialEquations(ListOfEquations&& eqns) :
    equations(std::forward<ListOfEquations>(eqns)),
    variables(linear_impl::get_variables<Variable>(equations)),
    sol(DefaultLinAlgAlgorithms::lsolve(lst{},variables))
=======
public:
  template< typename ListOfEquations>
  AbstractPolynomialEquations(ListOfEquations&& eqns) : 
    equations{expand(std::forward<ListOfEquations>(eqns))},
    variables{linear_impl::get_variables<Variable>(equations)}, 
    sol{DefaultLinAlgAlgorithms::lsolve(lst{},variables)} 
>>>>>>> e1b35850f1511f73a83d575c4b36cfdbbb0af0e4:linearsolve.h
  {
  }
  bool eliminate_linear_equations() {
		lst linear=linear_equations();
		if (linear==lst{}) return false;
    update_solution(DefaultLinAlgAlgorithms::lsolve(linear,variables));
    return true;
  }
  lst solution() const {return sol;}
  void dbgprint(ostream& os) const {
	  os<<"equations: "<<equations<<endl;
	  os<<"variables: "<<variables<<endl;
	  os<<"sol: "<<sol<<endl;
  }
};


template<class Variable,class Parameter>
class LinearEquationsWithParameters : protected AbstractPolynomialEquations<Variable> {
 	lst linear_equations() const override {	
		list<ex> linear;
		for (auto eq: this->equations) 
		  if (eq.is_polynomial(this->variables) && Degree<Poly<Variable>>(eq)<=1 && Degree<Poly<Parameter>>(eq)==0  && !eq.is_zero()) linear.push_back(eq==0);
		return linear;
	}
public:
	using AbstractPolynomialEquations<Variable>::AbstractPolynomialEquations;
	using AbstractPolynomialEquations<Variable>::eliminate_linear_equations;
	using AbstractPolynomialEquations<Variable>::solution;
	using AbstractPolynomialEquations<Variable>::dbgprint;
	lst always_solution() const {
		lst solution;
		list<ex> surviving_variables;
		for (auto eq : this->equations)
			GetSymbols<Variable>(surviving_variables,eq);
		lst to_zero;
		for (auto x: surviving_variables) to_zero.append(x==0);
		for (auto x: this->sol)
			solution.append(x.lhs()==x.rhs().subs(to_zero));
		for (auto x: to_zero)
			solution.append(x);
		return solution;
	}
};


template<class Variable>
class PolynomialEquations : protected AbstractPolynomialEquations<Variable>  {
 	lst linear_equations() const override {	
		list<ex> linear;
		for (auto eq: this->equations) 
		  if (eq.is_polynomial(this->variables) && Degree<Poly<Variable>>(eq)<=1 && !eq.is_zero()) linear.push_back(eq==0);		
		return linear;
	}
public:
	using AbstractPolynomialEquations<Variable>::AbstractPolynomialEquations;
	using AbstractPolynomialEquations<Variable>::eliminate_linear_equations;
	using AbstractPolynomialEquations<Variable>::solution;
	using AbstractPolynomialEquations<Variable>::dbgprint;
};


template<class Variable,class Parameter>
class LinearEquationsWithParameters : protected AbstractPolynomialEquations<Variable> {
 	lst linear_equations() const override {	
		list<ex> linear;
		for (auto eq: this->equations) 
		  if (eq.is_polynomial(this->variables) && Degree<Poly<Variable>>(eq)<=1 && Degree<Poly<Parameter>>(eq)==0  && !eq.is_zero()) linear.push_back(eq==0);
		return linear;
	}
public:
	using AbstractPolynomialEquations<Variable>::AbstractPolynomialEquations;
	using AbstractPolynomialEquations<Variable>::eliminate_linear_equations;
	using AbstractPolynomialEquations<Variable>::solution;
	lst always_solution() const {
		lst solution;
		list<ex> surviving_variables;
		for (auto eq : this->equations)
			GetSymbols<Variable>(surviving_variables,eq);
		lst to_zero;
		for (auto x: surviving_variables) to_zero.append(x==0);
		for (auto x: this->sol)
			solution.append(x.lhs()==x.rhs().subs(to_zero));
		for (auto x: to_zero)
			solution.append(x);
		return solution;
	}
};


template<class Variable>
class PolynomialEquations : protected AbstractPolynomialEquations<Variable>  {
 	lst linear_equations() const override {	
		list<ex> linear;
		for (auto eq: this->equations) 
		  if (eq.is_polynomial(this->variables) && Degree<Poly<Variable>>(eq)<=1 && !eq.is_zero()) linear.push_back(eq==0);		
		return linear;
	}
public:
	using AbstractPolynomialEquations<Variable>::AbstractPolynomialEquations;
	using AbstractPolynomialEquations<Variable>::eliminate_linear_equations;
	using AbstractPolynomialEquations<Variable>::solution;
};
}


template<typename Variable, typename ParametrizedClass, typename ListOfEquations, typename IsAdmissibleSolution> 
bool impose_polynomial_eqns(ParametrizedClass& parametrized_object, ListOfEquations&& eqns,const IsAdmissibleSolution& isAdmissibleSolution)
{
  linear_impl::PolynomialEquations<Variable> equations(std::forward<ListOfEquations>(eqns));
  while (equations.eliminate_linear_equations()) 
    if (equations.solution()==lst{}) return false;
  auto solution=equations.solution();
  nice_log<<solution<<endl;
  if (!isAdmissibleSolution(solution)) return false;
	exvector to_zero;
	std::transform(solution.begin(),solution.end(),std::back_insert_iterator<exvector>(to_zero),[](ex equation) {return equation.lhs()-equation.rhs();});
	parametrized_object.DeclareZero(to_zero.begin(),to_zero.end());
  return true;
}


}
#endif
