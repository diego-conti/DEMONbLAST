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
    
#ifndef ARROW_H
#define ARROW_H
#include "includes.h"

//represents an unlabeled arrow (in) -> (out)
struct Arrow {
	int node_in, node_out;
	void apply_permutation(const vector<int>& sigma) {node_in=sigma[node_in]; node_out=sigma[node_out];}
	bool operator==(const Arrow& other) const {return node_in==other.node_in && node_out==other.node_out;}
	void invert(int number_of_nodes) {node_in=number_of_nodes-node_in-1; node_out=number_of_nodes-node_out-1;}
};


//represents a labeled arrow (in) ->{^label} (out)
struct LabeledArrow {
	int node_in, node_out;
	int label = NO_LABEL;
	LabeledArrow() = default;
	LabeledArrow(const Arrow& arrow) {node_in=arrow.node_in; node_out=arrow.node_out;}

	void apply_permutation(const vector<int>& sigma) {
		node_in=sigma[node_in]; 
		node_out=sigma[node_out];
		if (has_label()) label=sigma[label];
	}
	bool operator==(const LabeledArrow& other) const {return node_in==other.node_in && node_out==other.node_out && label==other.label;}
	bool has_label() const {return label!=NO_LABEL;}
	void invert(int number_of_nodes) {
	  node_in=number_of_nodes-node_in-1; node_out=number_of_nodes-node_out-1;
	  if (has_label()) label=  number_of_nodes-label-1;
	}
private:
	static const int NO_LABEL=-1;
};


struct LabeledDoubleArrow {
  LabeledDoubleArrow(int from, int to, pair<int,int> label) : node_in{from}, node_out{to}, label{label} {}
	int node_in, node_out;
	pair<int,int> label;
	void apply_permutation(const vector<int>& sigma) {		
	  node_in=sigma[node_in]; 
		node_out=sigma[node_out];
		if (has_label()) label={sigma[label.first],sigma[label.second]};
  }
	bool operator==(const Arrow& other) const {return node_in==other.node_in && node_out==other.node_out;}
	bool has_label() const {return label!=pair<int,int>{};}
};

ostream& operator<<(ostream& os, Arrow arrow);
ostream& operator<<(ostream& os, LabeledArrow arrow);
ostream& operator<<(ostream& os, LabeledDoubleArrow arrow);

#endif
