#include "../weightbasis.cpp"
#include "../niceliegroup.cpp"
#include "../labeled_tree.cpp"
#include "../tree.cpp"
#include "../partitions.cpp"
#include "../gauss.cpp"
#include "../liegroupsfromdiagram.cpp"
#include "../filter.cpp"
#include "../diagramprocessor.h"
#include "../niceeinsteinliegroup.cpp"
#include "../permutations.cpp"
#include "../linearinequalities.h"

#include "dump.h"

void test_tree(LabeledTree tree) {
	tree.invert_nodes();
	CoefficientConfigurationWithoutRedundantParameter configuration{tree.weight_basis()};		
	stringstream os;
	while (configuration) {
		for (auto weight : configuration.weights())
				os<<weight.value<<weight<<endl;
		++configuration;
		os<<endl;
	}
	dump("signs"+tree.name(),os.str());
}

void test_tree(vector<int> partition, int hash) {
  auto diagrams = nice_diagrams(partition,Filter{});
	for (auto& diagram : diagrams) 
	  if (diagram.hash()==hash) {
	  	test_tree(diagram); return;
	  }
 	cerr<<"ERROR: no diagram in "<<horizontal(partition)<<" with hash "<<hash<<" found"<<endl;
}

void test_signs(int n) {
  stringstream os;
  for (auto conf:SignConfiguration::all_configurations(n))
  	os<<conf<<endl;
  dump("signs"+to_string(n),os.str());  
}

void test_inequalities() {
	StructureConstant a{N.a}, b{N.b};
	lst eqns;
	eqns=2*a+3*b+1, -a+b;
	LinearInequalities l1{eqns.begin(),eqns.end(),StructureConstant{}};
	assert((move(l1).has_solution()));
	eqns.remove_all();
	eqns = 2*a+2*b-1,a-b,b;
	LinearInequalities l2 {eqns.begin(),eqns.end(),StructureConstant{}};
	assert((move(l2).has_solution()));
	eqns.remove_all();
	eqns = 2*a+2*b-1,-a-b;
	LinearInequalities l3 {eqns.begin(),eqns.end(),StructureConstant{}};	
	assert((!move(l3).has_solution()));
}

int main() {
	test_inequalities();
	test_signs(0);
	test_signs(1);
	test_signs(2);
	test_signs(3);
	test_tree({4,2,2},3004);
	test_tree({4,3},82);
	test_tree({2,1,1,1,1,1,1},484836536);
}					
					
