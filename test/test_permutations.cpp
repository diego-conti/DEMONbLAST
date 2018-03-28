#include "../permutations.h"
#include <cassert>
#include <iostream>

void test_permutations(const vector<int>& from, const vector<int>& to) {
  for (auto sigma : PermutationsPreservingHash{from,to}.all_permutations())  {
    assert(sigma.size()==from.size());
    for (int i=0;i<sigma.size();++i)
      assert (to[sigma[i]]==from[i]);  
    for (int i=0;i<sigma.size();++i)
    for (int j=i+1;j<sigma.size();++j)
      assert(sigma[i]!=sigma[j]);
  }
}

void test_permutations() {
  test_permutations({1,2,3},{1,2,3});
  test_permutations({1,1,1},{1,1,1});
  test_permutations({1,2,1},{2,1,1});
  test_permutations({1,2,3,2,3},{3,3,2,2,1});
}

int main() {
  cout<<"testing permutations...";
  test_permutations();
  cout<<"OK"<<endl;
}
