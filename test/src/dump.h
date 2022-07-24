#include <fstream>

using namespace std;

template<typename Container>
void dump_container(string filename, Container&& trees) {
  ofstream file(filename+".out",std::ofstream::out | std::ofstream::trunc);	
  for (auto tree : trees)
   file<<tree;
}

template<typename Object>
void dump(string filename, const list<Object>& container) {
  dump_container(filename,container);
}
template<typename Object>
void dump(string filename, const vector<Object>& container) {
  dump_container(filename,container);
}

template<typename Object>
void dump(string filename, const Object& object) {
  ofstream file(filename+".out",std::ofstream::out | std::ofstream::trunc);	
  file<<object;
}



