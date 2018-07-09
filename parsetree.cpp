#include "parsetree.h"

namespace DiagramParser {

Token parse_token(istream& s) {
	if (s.peek()=='\n') {s.get(); return {TokenType::END};}
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
