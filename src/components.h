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
#ifndef COMPONENTS_H
#define COMPONENTS_H

#include "labeled_tree.h"

class Components {
	vector<set<int>> components;
	auto component_with(int i) {
		auto I=components.begin();
		while (I->find(i)==I->end()) ++I;
		return I;
	}
	void connect(int i, int j) {
		auto I=component_with(i), J=component_with(j);
		if (I!=J) {
			I->insert(J->begin(),J->end());
			components.erase(J);		
		}
	}
public:
	Components(const LabeledTree& diagram) {
		for (int i=0;i<diagram.number_of_nodes();++i) components.push_back({i});
		for (auto& weight : diagram.weights()) {		
			connect(weight.node_in1, weight.node_in2);
			connect(weight.node_in1,weight.node_out);
		}
	}
	int size() const {return components.size();}
};

#endif
