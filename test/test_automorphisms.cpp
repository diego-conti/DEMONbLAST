#include "../weightbasis.cpp"
#include "../niceliegroup.cpp"
#include "../labeled_tree.cpp"
#include "../tree.cpp"
#include "../partitions.cpp"
#include "../gauss.cpp"
#include "../liegroupsfromdiagram.cpp"
#include "../filter.cpp"
#include "../weightmatrix.cpp"
#include "../antidiagonal.cpp"
#include "../implicitmetric.cpp"
#include "../adinvariantobstruction.cpp"
#include "../parsetree.cpp"
#include "../automorphisms.cpp"
#include "dump.h"

using namespace std;

template<typename Tree>
string automorphisms_as_string(const Tree& tree) {
	list<string> result;
	for (auto sigma : tree.nontrivial_automorphisms())
		result.push_back("("+horizontal(incremented(sigma))+")");
	return horizontal(result);
}

void test_automorphisms(const vector<int>& partition) {
	list<string> output;
	for (auto& tree : nice_diagrams(partition,Filter{},DiagramDataOptions{}))
		output.push_back(tree.as_string()+"\n"+automorphisms_as_string(tree)+"\n");
	dump("automorphisms"+horizontal(partition,""),output);
}

void test_automorphisms() {
  test_automorphisms({6,1,1});
  test_automorphisms({5,2,1});
}

int main() {
  cout<<"testing automorphisms...";
  test_automorphisms();
  cout<<"OK"<<endl;
}
