#include "../tree.cpp"
#include "../log.cpp"
#include "../partitions.cpp"
#include <iostream>

auto p=partitions(8);

void speed_1() {
  for (auto partition : p) {
    auto result= incomplete_trees(partition,8);
    for (auto x: result) cout<<x<<endl;
 }
}

void speed_2() {
  for (auto partition : p) {
    auto result= complete_trees(partition,8);
    for (auto x: result) cout<<x<<endl;
  }
}

int main() {
  speed_1();
}
