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
#ifndef EXPRESSION_PARSER_H
#define EXPRESSION_PARSER_H

template<typename Parameter> class ExpressionParser {
	ex parameters;
	static string next_piece(const string& s, const string& separator, size_t& pos) {
		size_t next=s.find(separator,pos);
		if (next==string::npos) {
			auto result=s.substr(pos);
			pos=string::npos;
			return result;
		}
		else {
			auto result=s.substr(pos,next-pos);
			pos=next+separator.size();
			return result;
		}
	}
public:
	ExpressionParser() {
		lst parameters;
		for (auto v: NameRange(N.a,1,10))
			parameters.append(Parameter{v});
		for (auto v: NameRange(N.A,1,10))
			parameters.append(Parameter{v});
		for (ex x : {
			Parameter{N.a},Parameter{N.b},Parameter{N.c},Parameter{N.d},Parameter{N.e},Parameter{N.f},Parameter{N.g},Parameter{N.h},Parameter{N.i},
			Parameter{N.j},Parameter{N.k},Parameter{N.l},Parameter{N.m},Parameter{N.n},Parameter{N.o},Parameter{N.p},Parameter{N.q},Parameter{N.r},
			Parameter{N.s},Parameter{N.t},Parameter{N.u},Parameter{N.v},Parameter{N.w},Parameter{N.x},Parameter{N.y},Parameter{N.z},
			Parameter{N.A},Parameter{N.B},Parameter{N.C},Parameter{N.D},Parameter{N.E},Parameter{N.F},Parameter{N.G},Parameter{N.H},Parameter{N.I},
			Parameter{N.alpha},Parameter{N.beta},Parameter{N.gamma},Parameter{N.delta},Parameter{N.lambda},Parameter{N.mu},Parameter{N.xi},Parameter{N.chi}
		})
			parameters.append(x);
		this->parameters=parameters;
	}
	ex parse(const string& s) const {
		return ex{s,parameters};
	}
	exvector parse_vector(const string& s, const string& separator) const {
		if (s.empty()) return {};
		exvector result;
		size_t i=0;
		do {
			string piece=next_piece(s,separator,i);
			result.push_back(parse(piece));
		} while (i!=string::npos);
		return result;
	}
	exvector parse_vector_within_brackets(const string& s, const string& separator) const {
		auto i=s.find("["), j=s.rfind("]");
		if (i==string::npos || j==string::npos) throw std::invalid_argument("ExpressionParser: "+s+" is not a vector within []");		
		return parse_vector(s.substr(i+1,j-i-1),separator);
	}
};

#endif
