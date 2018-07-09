#include "../linearsolve.h"
#include <cassert>
#include <iostream>
#include <wedge/liegroup.h>

using namespace GiNaC;
using namespace std;
using namespace Wedge;

class TestLieGroup : public LieGroupHasParameters<true> , ConcreteManifold {
  StructureConstant a{N.a}, b{N.b}, c{N.c},d{N.d};
public:
  TestLieGroup() : ConcreteManifold{3} {
    Has_dTable::Declare_d(e(1),e(1)*e(2)+e(1)*e(3));
    Has_dTable::Declare_d(e(2),a*e(1)*e(3)+b*e(2)*e(3));  
    Has_dTable::Declare_d(e(3),c*e(1)*e(3)+d*e(1)*e(2));  
  }
};

void test7() {
  cout<<"testing PolynomialEquations...";
  lst eq;  
  StructureConstant a1(N.a(1)),a2(N.a(2)),a3(N.a(3)),a4(N.a(4)),a5(N.a(5)),a6(N.a(6)),a7(N.a(7)),a8(N.a(8));
  eq=a6-a4*a1,a5-a3+a7,a4-a2+a7,a5-a2,-1+a3,-1-a1,-1+a6;  
  linear_impl::PolynomialEquations<StructureConstant> equations{move(eq)};
  assert(equations.eliminate_linear_equations());
  assert(equations.eliminate_linear_equations());
  assert(!equations.eliminate_linear_equations());
  assert(ex{a5}.subs(equations.solution())==0);
  cout<<"OK"<<endl;
}

void test3() {
  cout<<"testing PolynomialEquations...";
  lst eq;  
  StructureConstant a(N.a), b(N.b), c(N.c);
  eq=a*b+a*c+a, b+c-1;
  assert(eq.op(0).is_polynomial(lst{a,b,c}));
  linear_impl::PolynomialEquations<StructureConstant> equations{move(eq)};
  assert(equations.eliminate_linear_equations());
  assert(equations.eliminate_linear_equations());
  assert(!equations.eliminate_linear_equations());
  assert(ex{a}.subs(equations.solution())==0);
  assert(ex{b+c-1}.subs(equations.solution())==0);
  cout<<"OK"<<endl;
}

int main() {
  test3();
   test7();  
  cout<<"testing impose_polynomial_eqns...";
  TestLieGroup G;
  lst eqns;
  G.GetEquations_ddZero(eqns);
  impose_polynomial_eqns<StructureConstant>(G,move(eqns),[](const lst&) {return true;});  
  G.Check_ddZero();
  cout<<"OK"<<endl;
}
