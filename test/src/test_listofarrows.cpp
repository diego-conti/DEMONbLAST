#include "tree.cpp"
#include "partitions.cpp"
#include <iostream>
#include "dump.h"

void test_ways() {
  dump("ways3",WaysToAddANode{3});
  dump("ways5",WaysToAddANode{5});
  dump("ways6",WaysToAddANode{6});
}

void test_iterator() {
  ListOfArrows list{2,Position::begin_t{}};
  cout<<"list 0"<<endl;
  for(auto i=list.begin(); i!=list.end();++i)
    i.print(cout); 
  for (auto x: list) cout<<x<<endl;
  ++list;
   cout<<"list 1"<<endl;
  for (auto x: list) cout<<x<<endl;
  ++list;
  cout<<"list 2"<<endl;
  for (auto x: list) cout<<x<<endl;
  ++list;
  cout<<"list 3"<<endl;
  for (auto x: list) cout<<x<<endl;
}
  
void test_list_of_arrows() {
  ListOfArrows list{2,Position::begin_t{}};
  assert(list.size()==2);
  assert(list.at(0)==false);
  assert(list.at(1)==false);
  assert(!list.is_this_the_end());
  ++list;
  assert(list.at(0)==true);
  assert(list.at(1)==false);
  assert(!list.is_this_the_end());
  ++list;
  assert(list.at(0)==false);
  assert(list.at(1)==true);
  assert(!list.is_this_the_end());
  ++list;
  assert(list.at(0)==true);
  assert(list.at(1)==true);
  assert(!list.is_this_the_end());
   ++list;
  assert(list.is_this_the_end());
}


int main() {
  test_list_of_arrows();
  test_ways();
}
