#include "../partitions.cpp"
#include <cassert>
#include <iostream>

bool is_one_of(vector<int> partition, list<vector<int>>& partitions) {
  return find(partitions.begin(),partitions.end(),partition)!=partitions.end();
}

void test_sorted_partitions() {
  auto part5=sorted_partitions(5);
  assert(part5.size()==7);
  assert(is_one_of({1,1,1,1,1},part5));
  assert(is_one_of({1,1,1,2},part5));
  assert(is_one_of({1,2,2},part5));
  assert(is_one_of({1,1,3},part5));
  assert(is_one_of({2,3},part5));
  assert(is_one_of({1,4},part5));
  assert(is_one_of({5},part5));
 }

void test_partitions() {
  auto part4=partitions(4);
  assert(part4.size()==8);
  assert(is_one_of({1,1,1,1},part4));
  assert(is_one_of({1,1,2},part4));
  assert(is_one_of({2,1,1},part4));
  assert(is_one_of({1,2,1},part4));
  assert(is_one_of({2,2},part4));
  assert(is_one_of({3,1},part4));
  assert(is_one_of({1,3},part4));
  assert(is_one_of({4},part4));  
 }


int main() {
  cout<<"testing partitions...";
  test_partitions();
  cout<<"OK"<<endl;
}
