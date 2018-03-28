#ifndef XGINAC_H
#define XGINAC_H

#include <ginac/ginac.h>

//allow range-based for loops for objects of type ex (or convertible to ex, like matrix and lst).

namespace GiNaC {
class const_ex_iterator {
	ex p;
	size_t index;
public:
	const_ex_iterator() = default;
	explicit const_ex_iterator(ex container, size_t idx=0) : p{container}, index{idx} {}
	ex operator*() const {return p.op(index);}
	const ex* operator->() const {return & (const_cast<ex&>(p).let_op(index));}
	const_ex_iterator& operator++() {++index; return *this;}
	const_ex_iterator operator++(int) {auto result=*this; ++index; return result;}
	const_ex_iterator& operator--() {--index; return  *this;}
	const_ex_iterator operator--(int) {auto result=*this; --index; return result;}
	bool operator!=(const const_ex_iterator& other) const {return p!=other.p || index!=other.index;}
	bool operator==(const const_ex_iterator& other) const {return p==other.p && index==other.index;}
	const_ex_iterator operator+=(int n) {index+=n; return *this;}
	const_ex_iterator operator-=(int n) {index-=n; return *this;}
	ex operator[] (int n) const {return p.op(index+n);}
	int operator-(const_ex_iterator j) const {return index-j.index;}
};

inline const_ex_iterator operator+ (const_ex_iterator i, int offset) {
	return i+=offset;
}
inline const_ex_iterator operator+ (int offset,const_ex_iterator i) {
	return i+=offset;
}
inline const_ex_iterator operator- (const_ex_iterator i, int offset) {
	return i-=offset;
}

inline const_ex_iterator begin(ex container) {return const_ex_iterator{container};}
inline const_ex_iterator end(ex container) {return const_ex_iterator{container,container.nops()};}

}

namespace std {
template<>
class iterator_traits<GiNaC::const_ex_iterator> {
public:
	using difference_type=int;
	using value_type = GiNaC::ex;
	using pointer= GiNaC::ex*;
	using reference = GiNaC::ex&;
	using iterator_category = random_access_iterator_tag;
};
}

#endif
