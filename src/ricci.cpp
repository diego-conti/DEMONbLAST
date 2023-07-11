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
#include "ricci.h"

ex ricci_tensor(const Manifold& G, matrix metric_on_coframe) {
	auto P=PseudoRiemannianStructureByMatrix::FromMatrixOnCoframe (&G,G.e(),metric_on_coframe);
	PseudoLeviCivitaConnection omega(&G,P);
	matrix ricci_tensor=omega.RicciAsMatrix();
	assert(ricci_tensor.cols()==G.Dimension());
	assert(ricci_tensor.rows()==G.Dimension());
	return ricci_tensor;
}

matrix diagonal_matrix(const exvector& diagonal) {
  matrix diagonal_matrix(diagonal.size(),diagonal.size());
  for (int i=0;i<diagonal.size();++i) diagonal_matrix(i,i)=diagonal[i];
  return diagonal_matrix;
}

ex ricci_tensor(const Manifold& G, exvector diagonal_metric) {
  assert(G.Dimension()==diagonal_metric.size());
	return ricci_tensor(G,diagonal_matrix(diagonal_metric));
}

ex ricci_operator(const Manifold& G, matrix metric_on_coframe) {
	return (metric_on_coframe* ricci_tensor (G,metric_on_coframe)).evalm();
}

ex ricci_operator(const Manifold& G, exvector diagonal_metric) {
	return ricci_operator(G,diagonal_matrix(diagonal_metric));
}
