#include "tree.cpp"
#include "partitions.cpp"
#include "labeled_tree.cpp"
#include "weightbasis.cpp"
#include "gauss.cpp"
#include "weightmatrix.cpp"
#include "antidiagonal.cpp"
#include "implicitmetric.cpp"
#include "adinvariantobstruction.cpp"
#include "parsetree.cpp"
#include "automorphisms.cpp"
#include "dump.h"


int main() {
    dump("incomplete21",incomplete_trees({2,1},3));
    dump("complete221",trees({2,2,1}));
    dump("complete21111",trees({2,1,1,1,1}));
    dump("complete321",trees({3,2,1}));
}
/*

void test_partitions() {
	for (auto partition : sorted_partitions(5))
		cout<<get_label(partition)<<endl;
}


void test_enlarge_tree() {
	Tree tree{3};
	auto trees=enlarge_tree(tree,2);
	for (auto t : trees) cout<<t;

}


struct TestLabeledTree {
	LabeledTree tree;
	TestLabeledTree(const Tree& t) :tree {t} {
		int node=tree.node_with_less_unlabeled_incoming_arrows();
		cout<<tree.node_with_less_unlabeled_incoming_arrows()<<endl;
		auto unlabeled_arrows = tree.unlabeled_arrows_pointing_to(node);
		for (auto x: unlabeled_arrows) cout<<x<<" ";
		auto with_label=tree.add_label(node,unlabeled_arrows[0],unlabeled_arrows[1]);
		if (unlabeled_arrows.size()>2)	with_label=with_label.add_label(node,unlabeled_arrows[0],unlabeled_arrows[2]);
		cout<<with_label;
		cout<<boolalpha<<with_label.satisfies_invariants_at_node(unlabeled_arrows[0])<< " "<<with_label.satisfies_invariants_at_node(unlabeled_arrows[1])<<endl;

		cout<<"labelings"<<endl;
		for (auto x: tree.labelings()) cout<<x;
	}
};

void test_complete_tree() {
	Tree tree{1};
	auto trees=complete_tree(tree,4);
	for (auto t : trees) {
		cout<<"tree"<<endl<<t;
		TestLabeledTree{t};
		cout<<endl;
	}

}


void test_hash() {
  for (auto tree : labeled_trees({4,2,1})) {
  	for (auto w : tree.weights()) cout<<w<<endl; {
      for (int i=0;i<tree.number_of_nodes();++i)
        cout<<"hash "<< i+1<< " = " <<std::hex <<tree.node_hash()[i] <<std::dec<<endl;
    }
    cout<<tree.hash()<<endl;
    }
}
*/
