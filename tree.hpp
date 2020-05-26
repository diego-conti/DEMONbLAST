
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
        os<<no_outgoing_concatenated_arrows[i]<< " outgoing concatenated arrows of length "<<i<<endl;
      for (int i=0;i<no_incoming_concatenated_arrows.size();++i)
        os<<no_incoming_concatenated_arrows[i]<< " incoming concatenated arrows of length "<<i<<endl;
    }
  protected:
    vector<int> no_outgoing_concatenated_arrows;	//TODO reserve size?
    vector<int> no_incoming_concatenated_arrows;
    static void resize_if_needed(vector<int>& v, int new_size) {
      if (v.size()<new_size) v.resize(new_size);
    }
	};

	class NodePathDataOrdered : public NodePathData {
	public:
    void add_incoming_arrows(const NodePathDataOrdered& incoming_node_data) {
      const vector<int>&  long_incoming_arrows= incoming_node_data.no_incoming_concatenated_arrows;
      add_long_arrows(no_incoming_concatenated_arrows,long_incoming_arrows);
    }
    void add_outgoing_arrows(const NodePathDataOrdered& outgoing_node_data) {
      const vector<int>& long_outgoing_arrows= outgoing_node_data.no_outgoing_concatenated_arrows;
      add_long_arrows(no_outgoing_concatenated_arrows,long_outgoing_arrows);
    }
    static void add_long_arrows(vector<int>& no_concatenated_arrows, const vector<int>& no_long_concatenated_arrows) {
      resize_if_needed(no_concatenated_arrows,no_long_concatenated_arrows.size()+1);
      ++no_concatenated_arrows[0];
      for (int i=0;i<no_long_concatenated_arrows.size();++i)
        no_concatenated_arrows[i+1]+=no_long_concatenated_arrows[i];
    }
	};

	class NodePathDataUnordered : public NodePathData {
	public:
		void add_outgoing_power_arrow(int length) {
			resize_if_needed(no_outgoing_concatenated_arrows,length);
			++no_outgoing_concatenated_arrows[length-1];
		}
		void add_incoming_power_arrow(int length) {
			resize_if_needed(no_incoming_concatenated_arrows,length);
			++no_incoming_concatenated_arrows[length-1];
		}
	};
	
	class Paths {
		using Path=list<int>;
		list<Path> paths;
		static Path join(Path path, int node) {
			path.push_back(node);
			assert(path.front()!=path.back());
			return path;
		}
		static Path join(int node, Path path) {
			path.push_front(node);
			assert(path.front()!=path.back());
			return path;
		}
		static Path join(Path path1, const Path& path2) {
			path1.insert(path1.end(),path2.begin(),path2.end());
			assert(path1.front()!=path1.back());
			return path1;
		}
	
		template<typename ArrowType>
		void add_arrow(ArrowType arrow) {
			list<Path> extendable_forward, extendable_backwards;
			for (auto path: paths)
				if (path.front()==arrow.node_out) extendable_backwards.push_back(path);
				else if (path.back()==arrow.node_in) extendable_forward.push_back(path);
			for (auto& path: extendable_forward)
				paths.push_back(join(path,arrow.node_out));
			for (auto& path: extendable_backwards)
				paths.push_back(join(arrow.node_in,path));
			for (auto& path1 : extendable_forward)
				for (auto& path2: extendable_backwards)
					paths.push_back(join(path1,path2));
			paths.push_back({arrow.node_in,arrow.node_out});
		}
		int size;	
	public:
		template<typename Tree>
		Paths(const Tree& tree) : size(tree.number_of_nodes()) {
			for (auto arrow : tree.arrows()) add_arrow(arrow);
		}

		vector<NodePathDataUnordered> node_data() const {
			vector<NodePathDataUnordered> result{size};		
		  for (auto& path: paths) {
				result[path.front()].add_outgoing_power_arrow(path.size()-1);
				result[path.back()].add_incoming_power_arrow(path.size()-1);
		 	}
			 return result;
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


template<typename Arrow> void TreeBase<Arrow>::compute_hash_and_cache_result_unordered() const {
	auto node_data=tree_impl::Paths {*this}.node_data();
  tree_hash=0;
  for (int node=0;node<no_nodes;++node)      
    tree_hash+=(nodes_hash[node]=node_data[node].hash());
  if (tree_hash==HASH_NOT_COMPUTED) tree_hash=~HASH_NOT_COMPUTED;
}	

template<typename Arrow> void TreeBase<Arrow>::compute_hash_and_cache_result_ordered() const {
    for (auto& arrow: arrows()) assert(arrow.node_in>arrow.node_out);
  	vector<tree_impl::NodePathDataOrdered> data(no_nodes);
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
  
  
