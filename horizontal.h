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
#ifndef HORIZONTAL_H
#define HORIZONTAL_H
#include <sstream>
#include <ginac/ex.h>
#include <ginac/operators.h>

using namespace std;
using namespace GiNaC;

inline string to_string_hex(unsigned int n) {
	stringstream ss;
  ss<<hex<<n;
  return ss.str();
}

template<typename Vector> string horizontal(const Vector& v, string sep=",") {
	if (v.empty()) return {};
	stringstream ss;
	auto i=v.begin();
	ss<<latex<<*i;
	while (++i!=v.end()) ss<<sep<<*i;
	return ss.str();
}

inline string to_string(const exvector& v) {
  return "("+horizontal(v)+")\n";
}

#endif
