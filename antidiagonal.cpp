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
#include "weightbasis.h"
#include <wedge/ginaclinalg.h>
#include <sstream>
#include "liegroupsfromdiagram.h"
#include "horizontal.h"
#include "weightmatrix.h"
#include "antidiagonal.h"


bool conflict(Weight weight,const WeightMatrix& weight_matrix) {
	return any_of(weight_matrix.weight_begin(),weight_matrix.weight_end(),
		[weight] (Weight weight2) {
			return (weight.node_in1==weight2.node_in1 || weight.node_in1==weight2.node_in2 || weight.node_in2==weight2.node_in1 ||weight.node_in2==weight2.node_in2 ) && weight.node_out==weight2.node_out
				|| (weight.node_in1==weight2.node_in1 && weight.node_in2==weight2.node_in2) 
				|| (weight.node_in2==weight2.node_in1 && weight.node_in1==weight2.node_in2);
		});
}

bool sigma_defines_ricci_flat(const WeightMatrix& weight_matrix,const OrderTwoAutomorphism& sigma) {
	for (auto weight=weight_matrix.weight_begin();weight!=weight_matrix.weight_end();++weight)
		if (conflict(sigma.apply(*weight),weight_matrix)) return false;
	return true;
}

bool cycle_of_length_two(Weight weight1, Weight weight2) {
	return (weight1.node_in1==weight2.node_out ||	weight1.node_in2==weight2.node_out) &&
		(weight2.node_in1==weight1.node_out ||	weight2.node_in2==weight1.node_out);
}

bool sigma_enhanced_has_cycle_of_length_two(const WeightMatrix& weight_matrix,const OrderTwoAutomorphism& sigma) {
	for (auto weight=weight_matrix.weight_begin();weight!=weight_matrix.weight_end();++weight)
	for (auto weight2=weight_matrix.weight_begin();weight2!=weight_matrix.weight_end();++weight2)
		if (cycle_of_length_two(*weight,sigma.apply(*weight2))) return true;
	return false;
}

list<OrderTwoAutomorphism> ricci_flat_sigma(const WeightMatrix& weight_matrix) {
	list<OrderTwoAutomorphism> result;
	int n=weight_matrix.cols();
	for (int k=1;2*k<=n;++k) 
		for (OrderTwoAutomorphism sigma{n,k}; sigma;++sigma) 
			if (sigma_defines_ricci_flat(weight_matrix,sigma))
				result.push_back(sigma);
	return result;
}

//
