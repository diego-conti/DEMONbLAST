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
#include "../weightmatrix.cpp"
#include "../antidiagonal.cpp"
#include "../implicitmetric.cpp"

#include "dump.h"

void test_table_mode(vector<int> partition,ostream& os) {
  auto processor = DiagramProcessor{lie_algebra_table};
  processor.invert_nodes();
  stringstream output;
  auto diagrams = nice_diagrams(partition,Filter{},DiagramDataOptions{});
  int count=0;
	for (auto& diagram: diagrams) {
      diagram.add_number_to_name(++count);
			os<<processor.process(diagram).data;  
	}
}

void test_table_mode(vector<int> partition) {
  stringstream os;
  test_table_mode(partition,os);
  dump("table"+horizontal(partition,""),os.str());  
}

int main() {
  test_table_mode({2,1,1,1});
  test_table_mode({2,1,1,1,1});
  test_table_mode({2,1,1,1,1,1});
  test_table_mode({2,1,1,1,1,1,1});
//  test_table_mode({2,1,1,1,1,1,1,1});
}					
					
