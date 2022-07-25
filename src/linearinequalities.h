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
#ifndef LINEARINEQUALITIES_H
#define LINEARINEQUALITIES_H

<<<<<<< HEAD:src/linearinequalities.h
#include "includes.h"
=======
#include <ginac/ginac.h>
>>>>>>> e1b35850f1511f73a83d575c4b36cfdbbb0af0e4:linearinequalities.h

using namespace GiNaC;

template<typename iter_begin,typename iter_end, typename Closure>
iter_begin minimize_value(iter_begin begin, iter_end end, Closure&& closure) {
	int minimum=numeric_limits<int>::max();
	auto minimum_iter=begin;
	while (begin!=end) {
		int value=closure(*begin);
		if (value<minimum) {
			minimum=value;
			minimum_iter=begin;
		}
		++begin;
	}
	return minimum_iter;
}

template<typename iter_begin,typename iter_end, typename T, typename Closure>
T transform_accumulate(iter_begin begin, iter_end end, T value, Closure&& closure) {
	while (begin!=end) value+=closure(*begin++);
	return value;
}

class LinearInequalities {
	list<ex> inequalities;
	list<ex> unknowns;
	int count_occurrences(ex variable) const {
		return transform_accumulate(inequalities.begin(),inequalities.end(),0,[variable] (ex inequality) {
			exset found;
			return inequality.find(variable,found)? 1 : 0;
		});
	}
	ex variable_that_appears_in_the_least_equations() const {
		return *minimize_value(unknowns.begin(),unknowns.end(),[this] (ex variable) {
			int occurrences=count_occurrences(variable);
			return occurrences? occurrences : numeric_limits<int>::max();
		});
	}
	
	bool remove_constant_inequalities() {
		auto iter=inequalities.begin();
		while (iter!=inequalities.end()) {
			if (is_a<numeric>(*iter)) {
				if (*iter>0) iter=inequalities.erase(iter);
				else return false;
			}
			else ++iter;
		}	
		return true;
	}
	
	void eliminate() {
		if (!unknowns.empty())
			eliminate(variable_that_appears_in_the_least_equations());
	}
	void eliminate(ex x) {
		list<ex> greater_than_x, smaller_than_x, not_depending_on_x;
		for (ex inequality: inequalities) {
			ex coeff=inequality.coeff(x);
			assert(is_a<numeric>(coeff));
			if (coeff.is_zero()) not_depending_on_x.push_back(inequality);
			else if (coeff>0) smaller_than_x.push_back(-inequality/coeff);
			else greater_than_x.push_back(-inequality/coeff);
		}
		inequalities=not_depending_on_x;
		for (auto greater : greater_than_x)
		for (auto smaller : smaller_than_x)
			inequalities.push_back((greater-smaller).expand());
	}
public:
	template<typename Iterator,typename Variable>
	LinearInequalities(Iterator begin, Iterator end,Variable) : inequalities(begin,end) {
		GetSymbols<Variable>(unknowns,begin,end);
	}

	//applies the Fourier-Motzkin algorithm to determine whether the inequalities have a common solution
	bool has_solution() && {
		if (!remove_constant_inequalities()) return false;
		while (!inequalities.empty()) {
			if (!remove_constant_inequalities()) return false;
			eliminate();
		}
		return true;
	}
};

#endif
