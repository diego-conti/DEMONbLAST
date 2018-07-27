#include "../weightbasis.cpp"
#include "../niceliegroup.cpp"
#include "../labeled_tree.cpp"
#include "../tree.cpp"
#include "../partitions.cpp"
#include "../gauss.cpp"
#include "../liegroupsfromdiagram.cpp"
#include "../filter.cpp"
#include "../weightmatrix.cpp"

#include "dump.h"

void test_weight_basis() {
	auto trees = nice_diagrams({3,2,1},Filter{});
	for (auto tree : trees) {
		cout<<tree;
		for (auto weight_and_coefficient : WeightBasis{tree}.weights_and_coefficients())
			cout<<weight_and_coefficient<<endl;
	}
}

auto lie_groups(vector<int> partition) {
  auto result = nice_diagrams(partition,Filter{});
		for (auto& tree : result) {
					string groups;
					for (auto group : NiceLieGroup::from_weight_basis(WeightBasis{tree})) 
						groups+=to_string(group)+"\n";						
					if (groups.empty()) {
						groups="no Lie algebra";
					}
					tree.set_name(groups);					
				}
		return result;
}

int main() {
  dump("groups212",lie_groups({2,1,2}));
  dump("groups41",lie_groups({4,1}));
  dump("groups3111",lie_groups({3,1,1,1}));
}					
					
