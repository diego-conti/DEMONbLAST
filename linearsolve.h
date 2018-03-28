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

#include <wedge/parameters.h>
#include <wedge/expressions.h>
#include <wedge/lambda.h>
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



template<class Variable>
class PolynomialEquations {
  lst equations;
  lst variables;
  lst sol;
  
  lst linear_equations() const {
		list<ex> linear;
		for (auto eq: equations) 
		  if (eq.is_polynomial(variables) && Degree<Poly<Variable>>(eq)<=1 && !eq.is_zero()) linear.push_back(eq==0);		
		return linear;
	}
  
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
    for (int i=0;i<eqns.nops();++i)
      eqns.let_op(i)=eqns.op(i).expand();
    return std::move(eqns);
  }
public:
  template< typename ListOfEquations>
  PolynomialEquations(ListOfEquations&& eqns) : 
    equations{expand(std::forward<ListOfEquations>(eqns))},
    variables{linear_impl::get_variables<Variable>(equations)}, 
    sol{DefaultLinAlgAlgorithms::lsolve(lst{},variables)} 
  {
  }
  bool eliminate_linear_equations() {
		lst linear=linear_equations();
		if (linear==lst{}) return false;
    update_solution(DefaultLinAlgAlgorithms::lsolve(linear,variables));
    return true;
  }
  lst solution() const {return sol;}
};
}


template<typename Variable, typename ParametrizedClass, typename ListOfEquations, typename IsAdmissibleSolution> 
bool impose_polynomial_eqns(ParametrizedClass& parametrized_object, ListOfEquations&& eqns,const IsAdmissibleSolution& isAdmissibleSolution)
{
  linear_impl::PolynomialEquations<Variable> equations{std::forward<ListOfEquations>(eqns)};
  while (equations.eliminate_linear_equations()) 
    if (equations.solution()==lst{}) return false;
  auto solution=equations.solution();
  nice_log<<solution<<endl;
  if (!isAdmissibleSolution(solution)) return false;
	exvector to_zero;
	std::transform(solution.begin(),solution.end(),back_insert_iterator<exvector>(to_zero),[](ex equation) {return equation.lhs()-equation.rhs();});
	parametrized_object.DeclareZero(to_zero.begin(),to_zero.end());
  return true;
}


}
#endif
