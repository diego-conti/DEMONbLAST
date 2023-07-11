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
#ifndef ADINVARIANTOBSTRUCTION_H
#define ADINVARIANTOBSTRUCTION_H
#include "weightbasis.h"
#include "labeled_tree.h"
#include "includes.h"

//return true if the obstruction of Prop 2.2 in  V. del Barco, Lie algebras admitting symmetric, invariant andnondegenerate bilinear forms, arxiv 1602.08286 is satisfied
bool passes_obstruction_for_ad_invariant_metric(const LabeledTree& diagram);

#endif
