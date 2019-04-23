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
	ParseError(const string& error) : runtime_error(error) {};
};

Token parse_token(istream& s) {
	if (s.peek()=='\n') return {TokenType::END};
	char ch = 0;
	int value;
	s>>ch;
	switch (ch) {
		case 0:
		return {TokenType::END};
		case '0':
		case '1':
		case '2':
		case '3':
		case '4':
		case '5':
		case '6':
		case '7':
		case '8':
		case '9':
			s.putback(ch);
			s>>value;
			return {TokenType::NUMBER,value};
		case '[':
			return {TokenType::OPEN_BRACKET};
		case ']':
			return {TokenType::CLOSE_BRACKET};
		case ':':
		case ',':
		case ';':
			return {TokenType::SEPARATOR};		
		case '-':
			s>>ch;
			if (ch=='>') return {TokenType::ARROW};
			else throw ParseError(string{"> expected after - instead of "}+ch);
		default:
			throw ParseError("Unexpected char: "+ch);
		}
}

int parse_number(istream& s) {
	auto token=parse_token(s);
	if (token.type!=TokenType::NUMBER) throw ParseError("Number expected");
	return token.data;
}

void skip_token(istream& s, TokenType type) {
	auto token=parse_token(s);
	if (token.type!=type) throw ParseError(Token{type}.to_string()+ " expected, found "+token.to_string());
}

bool parse_arrow(istream& s, LabeledArrow& arrow) {
		auto token= parse_token(s);
		if (token.type==TokenType::END) return false;
		if (token.type!=TokenType::SEPARATOR) throw ParseError("Separator expected, found "+token.to_string());
		arrow.node_in=parse_number(s)-1;
		skip_token(s,TokenType::ARROW);
		skip_token(s,TokenType::OPEN_BRACKET);
		arrow.label=parse_number(s)-1;				
		skip_token(s,TokenType::CLOSE_BRACKET);
		arrow.node_out=parse_number(s)-1;
		return true;
}

}

#endif
