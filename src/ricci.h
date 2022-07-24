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
#ifndef RICCI_H
#define RICCI_H

#include "includes.h"

using namespace Wedge;

class ScalarProductDefinedByMatrix : public BilinearFormWithFrame {
	matrix m_, m_inverse_;
public:
	ScalarProductDefinedByMatrix(const Frame& frame, matrix m) : BilinearFormWithFrame(frame),m_{m}, m_inverse_{m.inverse()} {}
protected:
	ex MatrixEntry(OneBased i, OneBased j) const {
		return m_(i-1,j-1);
	}
	ex InverseMatrixEntry(OneBased i, OneBased j) const {
		return m_inverse_(i-1,j-1);
	}
};


class SomePseudoRiemannianStructure : public PseudoRiemannianStructure {
	const ScalarProductDefinedByMatrix scalar_product;
public:
/** @brief 
 *  @param manifold The manifold on which the structure is defined.
 *  @param orthonormal_frame A coframe with respect to which the metric is defined; could be orthonormal or not.
 *  @param p An integer that determines the signature
 *
 * The orthonormal frame e_1,..., e_n is assumed to satisfy <e_i,e_j>=\delta_{ij} for i\leq p, -\delta_{ij} for i>p
*/
	SomePseudoRiemannianStructure(const Manifold* manifold, const Frame& orthonormal_frame, const matrix& m) : 
		GStructure{manifold, orthonormal_frame}, PseudoRiemannianStructure(manifold,orthonormal_frame), scalar_product{orthonormal_frame, m} {}
	const BilinearForm& ScalarProduct() const {return scalar_product;}

	pair<ex,matrix> DecomposeRicci(matrix ricci) const override {return {};}
};

ex ricci_tensor(const Manifold& G, matrix metric_on_coframe);

ex ricci_tensor(const Manifold& G, exvector diagonal_metric);

ex ricci_operator(const Manifold& G, matrix metric_on_coframe);

ex ricci_operator(const Manifold& G, exvector diagonal_metric);

#endif
