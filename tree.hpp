
template<typename ArrowType>
string TreeBase<ArrowType>::to_dot_string(string extra_data) const {
    stringstream os;
		os<<"digraph {"<<endl;
		os<<"label = \""+name()+extra_data+"\";"<<endl;
		for (int i=0;i<number_of_nodes();++i) os<<i+1<<endl;
		for (auto arrow: arrows()) os<<arrow;
		os<<"}"<<endl;
		return os.str();
}


namespace tree_impl {
	class NodePathData {
	public:
    void add_incoming_arrows(const NodePathData& incoming_node_data) {
      const vector<int>&  long_incoming_arrows= incoming_node_data.no_incoming_concatenated_arrows;
      add_long_arrows(no_incoming_concatenated_arrows,long_incoming_arrows);
    }
    void add_outgoing_arrows(const NodePathData& outgoing_node_data) {
      const vector<int>& long_outgoing_arrows= outgoing_node_data.no_outgoing_concatenated_arrows;
      add_long_arrows(no_outgoing_concatenated_arrows,long_outgoing_arrows);
    }
    int hash() const {    
      int result =no_outgoing_concatenated_arrows.size();
      for (int i=0;i<no_outgoing_concatenated_arrows.size();++i)
        result=result<<4 | no_outgoing_concatenated_arrows[i];
      for (int i=0;i<no_incoming_concatenated_arrows.size();++i)
        result=result<<4 | no_incoming_concatenated_arrows[i]/2;  //these are always even numbers so we save a bit by dividing by two
      return (result>0)? result : -result;
    }
    void print(ostream& os) const {
      for (int i=0;i<no_outgoing_concatenated_arrows.size();++i)
        cout<<no_outgoing_concatenated_arrows[i]<< " outgoing concatenated arrows of length "<<i<<endl;
      for (int i=0;i<no_incoming_concatenated_arrows.size();++i)
        cout<<no_incoming_concatenated_arrows[i]<< " incoming concatenated arrows of length "<<i<<endl;
    }
  private:
    vector<int> no_outgoing_concatenated_arrows;
    vector<int> no_incoming_concatenated_arrows;
    static void resize_if_needed(vector<int>& v, int new_size) {
      if (v.size()<new_size) v.resize(new_size);
    }
    static void add_long_arrows(vector<int>& no_concatenated_arrows, const vector<int>& no_long_concatenated_arrows) {
      resize_if_needed(no_concatenated_arrows,no_long_concatenated_arrows.size()+1);
      ++no_concatenated_arrows[0];
      for (int i=0;i<no_long_concatenated_arrows.size();++i)
        no_concatenated_arrows[i+1]+=no_long_concatenated_arrows[i];
    }
	};

}

template<typename Arrow> bool TreeBase<Arrow>::is_equivalent_to(const TreeBase& tree) const {
	if (hash()!=tree.hash()) return false;
	if (no_nodes!=tree.no_nodes) return false;
	if (arrows().size()!=tree.arrows().size()) return false;
	auto permutations =PermutationsPreservingHash{node_hash(), tree.node_hash()}.all_permutations(); 
	return any_of(permutations.begin(),permutations.end(),
	  [&tree,this] (const vector<int>& sigma) {return matches(tree,sigma);}
	  );
}


template<typename Arrow> list<vector<int>> TreeBase<Arrow>::nontrivial_automorphisms() const {
  list<vector<int>> automorphisms;
	auto permutations =PermutationsPreservingHash{node_hash(), node_hash()}.all_permutations(); 
  copy_if(permutations.begin(),permutations.end(),back_inserter(automorphisms),
    [this] (const vector<int>& sigma) {return matches(*this,sigma);}
  );
  automorphisms.remove(PermutationsPreservingHash::trivial_permutation(number_of_nodes()));
  return automorphisms;
}


template<typename Arrow> bool TreeBase<Arrow>::matches(const TreeBase& tree, const vector<int>& permutation) const {
	for (Arrow arrow: arrows()) {
		arrow.apply_permutation(permutation);
		if (find(tree.arrows().begin(),tree.arrows().end(),arrow)==tree.arrows().end()) return false;
	}
	return true;
}

template<typename Arrow> void TreeBase<Arrow>::compute_hash_and_cache_result() const {
  	vector<tree_impl::NodePathData> data(no_nodes);
    for (int node=no_nodes-1;node>=0;--node)
      for (auto& arrow: arrows()) 
         if (arrow.node_out==node) data[node].add_incoming_arrows(data[arrow.node_in]);
    for (int node=0;node<no_nodes;++node)
      for (auto& arrow: arrows()) 
         if (arrow.node_in==node) data[node].add_outgoing_arrows(data[arrow.node_out]);
    tree_hash=0;
    for (int node=0;node<no_nodes;++node)      
      tree_hash+=(nodes_hash[node]=data[node].hash());
    if (tree_hash==HASH_NOT_COMPUTED) tree_hash=~HASH_NOT_COMPUTED;
	}



template<typename Tree> void remove_equivalent_trees_same_hash(list<Tree>& list_of_trees) {
 for (auto i=list_of_trees.begin();i!=list_of_trees.end();++i) {
  auto j=i;
  ++j;
  while (j!=list_of_trees.end())
    if (i->is_equivalent_to(*j)) 
        j=list_of_trees.erase(j);
      else ++j;
 }
}


template<typename Tree> void remove_equivalent_trees(list<Tree>& list_of_trees) {
  map<int, list<Tree>> tree_by_hash;
  for (auto i=list_of_trees.begin();i!=list_of_trees.end();) {
    auto next=i; ++next;
    auto& list = tree_by_hash[i->hash()];
    list.splice(list.end(), list_of_trees, i);
    i=next;
   }
  assert(list_of_trees.empty());
  for (auto& hash_and_list : tree_by_hash) {
    remove_equivalent_trees_same_hash(hash_and_list.second);
    list_of_trees.splice(list_of_trees.end(), hash_and_list.second);
  }
  for (auto tree : list_of_trees)
    assert(tree.is_hash_computed());
}
  
  
