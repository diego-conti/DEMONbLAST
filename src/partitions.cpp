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
#include "partitions.h"
#include <sstream>
#include <set>
#include <algorithm>
#include <cassert>

using namespace std;

string get_label(vector<int> partition, string separator) {
  assert(partition.size()>0);
	stringstream s;
	auto i=partition.begin();
	s<<*i;
	while (++i!=partition.end()) s<<separator<<*i;
	return s.str();
}


list<vector<int>> sorted_partitions (int n, int maximum_width) {
	list<vector<int>> result;
	if (n==0) return {{}};
	for (int i=1;i<=n && i<=maximum_width;++i) {
		auto partitions_n_minus_one = sorted_partitions(n-i,i);
		for (auto& partition : partitions_n_minus_one) {
			partition.push_back(i);
			result.push_back(partition);
		}	
	}
	return result;
}

list<vector<int>> sorted_partitions (int n) {
	return sorted_partitions(n,n);
}

list<vector<int>> partitions (int n) {
	list<vector<int>> sorted=sorted_partitions(n);
	list<vector<int>> result;
	for (auto& partition : sorted) {
		set<vector<int>> permutations;
		vector<int> permuted = partition;
		do permutations.insert(permuted);			
		while (next_permutation(permuted.begin(),permuted.end()));		
		result.insert(result.end(),permutations.begin(),permutations.end());
	}
	return result;
}

