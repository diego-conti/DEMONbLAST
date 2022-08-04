#ifndef XGINAC_H
#define XGINAC_H

#include "includes.h"

inline GiNaC::exvector coefficients_of(GiNaC::exvector v, GiNaC::ex parameter) {
	transform(v.begin(),v.end(),v.begin(),[&parameter] (GiNaC::ex x) {return x.coeff(parameter);});
	return v;
}

#endif
