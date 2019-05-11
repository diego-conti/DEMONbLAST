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
#ifndef PERMUTATIONS_H
#define PERMUTATIONS_H

#include "includes.h"

class PermutationsPreservingHash {
  enum {ASSIGNED=-1};
  vector<int> hash_codes_from;
  vector<int> hash_codes_to;
  vector<int> current_permutation;
  
  void assign(list<vector<int>>& permutations, int first_unassigned, int to) {    
      int hash=hash_codes_to[to];
      if (hash_codes_from[first_unassigned]==hash) {
        hash_codes_to[to]=ASSIGNED;
        current_permutation[first_unassigned]=to;
        permutations.splice(permutations.end(),all_permutations(first_unassigned+1));        
        hash_codes_to[to]=hash;
      }  
  }
  list<vector<int>> all_permutations(int first_unassigned) {
    if (first_unassigned==hash_codes_from.size()) return {current_permutation};
    list<vector<int>> result;
    for (int assigned_to=0;assigned_to<hash_codes_to.size();++assigned_to) 
      assign(result, first_unassigned,assigned_to);
    return result;        
  }
public:
  PermutationsPreservingHash(const vector<int>& from, const vector<int>& to) : hash_codes_from(from), hash_codes_to(to),current_permutation(from.size()) {}
  list<vector<int>> all_permutations() {
    return all_permutations(0);
  }
  static vector<int> trivial_permutation(int dimension) {
    vector<int> result(dimension);
    iota(result.begin(),result.end(),0);
    return result;
  }
};

string permutation_to_string(vector<int> sigma);

#endif
