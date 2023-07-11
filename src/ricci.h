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
#ifndef RICCI_H
#define RICCI_H

#include "includes.h"

using namespace Wedge;


ex ricci_tensor(const Manifold& G, matrix metric_on_coframe);

ex ricci_tensor(const Manifold& G, exvector diagonal_metric);

ex ricci_operator(const Manifold& G, matrix metric_on_coframe);

ex ricci_operator(const Manifold& G, exvector diagonal_metric);

#endif
