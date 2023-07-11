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
#ifndef PARSE_TREE_H
#define PARSE_TREE_H

//5: 1->[2] 3, 1->[3] 4 means (0,0,12,13,0)
//5 means (0,0,0,0,0). Notice the absence of a colon.
#include "arrow.h"
#include <exception>

namespace DiagramParser {

enum class TokenType : char {
	NUMBER='0', ARROW='>', SEPARATOR=',', OPEN_BRACKET='[', CLOSE_BRACKET=']',END='E'
};


struct Token {
	TokenType type;
	int data;
	Token(TokenType type=TokenType::END, int data=-1) : type{type},data{data} {}
	string to_string() const {
		return (type==TokenType::NUMBER)? std::to_string(data) : string{static_cast<char>(type)};
	}
};

class ParseError : public std::runtime_error {
public:
	ParseError(const string& error) : std::runtime_error(error) {};
};

Token parse_token(istream& s);
bool parse_arrow(istream& s, LabeledArrow& arrow);

}

#endif
