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
#ifndef DEFORMATION_H
#define DEFORMATION_H
#include "linearsolve.h"
#include "weightbasis.h"


inline ex three_way_product(const exvector& x1, const exvector& x2, const exvector& x3) {
  ex result;
  assert(x1.size()==x2.size());
  assert(x1.size()==x3.size());
  for (int i=0;i<x1.size();++i)
    result+=x1[i]*x2[i]*x3[i];
  return result;
}

template<typename Parameter> 
exvector generic_vector(int length) {
  exvector result;
 	for (int i=1;i<=length;++i)
 	  result.push_back(Parameter{N.a(i)});
 	return result;
}

string lst_to_string(const lst& list) {
  stringstream s;
  s<<list;
  return s.str();
}

WEDGE_DECLARE_NAMED_ALGEBRAIC(DeformationParameter,realsymbol);

class Deformation : public HasParameters<DeformationParameter> {
  exvector M_Delta_v;
  lst deformation_equations;
  string data;
  
  string data_for_admissible_deformation() const {
        if (all_of(M_Delta_v.begin(),M_Delta_v.end(),[](ex x) {return x.is_zero();})) return "no nontrivial deformation";
        return horizontal(M_Delta_v) + " if "+ lst_to_string(deformation_equations);  
  }
  string data_for_no_admissible_deformation() const {
        return "no admissible deformation"; //al momento non succede mai perchè la deformazione zero è sempre ammissibile 
  }
public:
  Deformation(const WeightMatrix& weight_matrix)   {
    auto deformation_parameters = generic_vector<DeformationParameter>(weight_matrix.M_Delta().cols());
    M_Delta_v=weight_matrix.M_Delta().image_of(deformation_parameters);
    if (!weight_matrix.X_ijk().empty()) {
      for (int j=0;j<weight_matrix.M_Delta().cols(); ++j)   
        deformation_equations.append(three_way_product(weight_matrix.X_ijk(),M_Delta_v, weight_matrix.M_Delta().column(j)));  
      data = (impose_polynomial_eqns<DeformationParameter>(*this, deformation_equations,[](const lst&) {return true;}))?
        data_for_admissible_deformation() : data_for_no_admissible_deformation();
   }
  }    
  void DeclareConditions(const lst& list_of_equations) override {
      for (int i=0;i<deformation_equations.nops();++i)
        deformation_equations.let_op(i)=deformation_equations.op(i).subs(list_of_equations);   
      for (auto& param: M_Delta_v) param=param.subs(list_of_equations);   
  }
  string to_string() const {return data;}
};


#endif
