#ifndef DIAGRAMANALYZER_H
#define DIAGRAMANALYZER_H
#include <iostream>
#include <list>
#include <ginac/ginac.h>
#include "xginac.h"
#include <algorithm>
#include "labeled_tree.h"


//Represents a path in a diagram that does not strictly contain a cycle
//the diagram is assumed to have no loops (i.e. 1-cycles)
class Path {
	vector<int> sequence;
	vector<int> labels;
	Path enlarge(Weight weight) const {
		if (weight.node_in2==sequence.back()) swap(weight.node_in1, weight.node_in2);
		Path result=*this;
		result.sequence.push_back(weight.node_out);
		result.labels.push_back(weight.node_in2);
		return result;
	}
	bool can_be_joined_to(Weight weight) const {
		return sequence.empty() || sequence.back()==weight.node_in1 ||sequence.back()==weight.node_in2;
	}
	bool strictly_contains_a_cycle() const {		
		for (int i=1;i<sequence.size()-1;++i)
			if (sequence[i]==sequence.back()) return true;
		return false;
	}
	template<typename Iterator>
	list<Path> short_paths(Iterator weight_begin, Iterator weight_end) const {
		list<Path> result;
		for (;weight_begin!=weight_end;++weight_begin) {
			result.push_back(Path{weight_begin->node_in1,weight_begin->node_in2,weight_begin->node_out});
			result.push_back(Path{weight_begin->node_in2,weight_begin->node_in1,weight_begin->node_out});
		}		
		return result;
	}
	template<typename Iterator>
	list<Path> one_node_enlargements(Iterator weight_begin, Iterator weight_end) const {
		list<Path> result;
		for (;weight_begin!=weight_end;++weight_begin) {
			if (can_be_joined_to(*weight_begin)) {
				auto enlarged=enlarge(*weight_begin);
				if (!enlarged.strictly_contains_a_cycle()) result.push_back(enlarged);
			}
		}
		return result;
	}
	template<typename Iterator>
	list<Path> enlargements(Iterator weight_begin, Iterator weight_end) const {
		list<Path> result;
		auto enlarged= sequence.size()? one_node_enlargements(weight_begin,weight_end) : short_paths(weight_begin,weight_end);
		for (auto& path : enlarged) 
			result.splice(result.end(),path.enlargements(weight_begin,weight_end));		
		result.splice(result.end(),enlarged);
		return result;	
	}
	Path(int node_in, int label, int node_out) : sequence{node_in,node_out}, labels{label} {}
public:
	Path()=default;
	template<typename Iterator>
	static list<Path> paths(Iterator weight_begin, Iterator weight_end) {
		return Path{}.enlargements(weight_begin,weight_end);
	}
	int start_point() const {return sequence.front();}
	int end_point() const {return sequence.back();}
	string to_string() const {		
		return horizontal(incremented(sequence)," ");
	}
	bool cycle() const {return sequence.size()>1 && start_point()==end_point();}
	vector<int> cycle_labels_or_empty() const {return cycle()? labels : vector<int> {};}
};

inline ostream& operator<<(ostream& os,const Path& path) {
	return os<<path.to_string();
}

/*
class Cycle {
	Path path1, path2;
	bool touches(Weight weight) {
	
		}
public:
	Cycle(const Path& path1, const Path& path2) : path1{path1},path2{path2} {}
	
	template<typename Iterator>
	int extra_connected_weights(Iterator weight_begin, Iterator weight_end) {
		return count_if(weight_begin,weight_end,
			[this] (auto& weight) {return weight.node_out
	}
};
*/

class DiagramAnalyzer {
	list<Path> paths;	
	vector<int> nodes_for_lorentzian_diagonal_ricci_flat;
public:
	DiagramAnalyzer()=default;
	template<typename Weights>
	DiagramAnalyzer(int nodes, Weights&& weights) : paths{Path::paths(weights.begin(),weights.end())} {
		// Lorentzian diagonal ricci-flat condition	
		vector<int> result;
		auto centre=centre_basis(nodes, weights);
		auto gprimeperp=gprimeperp_basis(nodes, weights);
		for (int i=0;i<nodes;++i)
			if (all_of(centre.begin(),centre.end(),
				[&weights,i](int v_in_centre) {return reachable(weights,v_in_centre, i);})
			 	&& all_of(gprimeperp.begin(),gprimeperp.end(),
						[&weights,i](int v_in_gprimeperp) {return reachable(weights,i,v_in_gprimeperp) || !commutes(weights,i,v_in_gprimeperp);})
			)
				nodes_for_lorentzian_diagonal_ricci_flat.push_back(i);
	}
	
	/*
	list<Cycle> cycles() {
		list<pair<Path>> result;
		for (auto i=paths.begin();i!=paths.end();++i) {
			auto j=i;
			while (++j!=paths.end())
				if (i->start_point()==j->end_point() ||j->start_point()==i->end_point()) 
					result.emplace_back(*i,*j);
		}
	}*/
	list<Path> cycles() {
		list<Path> result;
		for (auto path : paths) 
			if (path.start_point()==path.end_point()) result.push_back(path);
		return result;
	}
	string killing() const {
		string result="label sequences in cycles:";
		set<vector<int>> labels;
		for (auto path : paths) 
			labels.insert(path.cycle_labels_or_empty());
		for (auto l: labels) 
			if (!l.empty()) result+=horizontal(incremented(l)," ")+"; ";
		return result;
	}	
	string analysis() const {
		return "Paths: "+horizontal(paths)+"\n"+"Lorentzian diagonal ricci-flat nodes: "+horizontal(incremented(nodes_for_lorentzian_diagonal_ricci_flat));
	}

	
private:
	template<typename Weights>
	static vector<int> centre_basis(int nodes, Weights&& weights) {
		auto basis=consecutive_numbers(0,nodes);
		for (auto weight : weights) {
			basis.erase(weight.node_in1);
			basis.erase(weight.node_in2);
		}
		return {basis.begin(),basis.end()};
	}

	template<typename Weights>
	static vector<int> gprimeperp_basis(int nodes, Weights&& weights) {
		auto basis=consecutive_numbers(0,nodes);
		for (auto weight : weights) basis.erase(weight.node_out);
		return {basis.begin(),basis.end()};
	}

	template<typename Weights>
	static bool reachable(Weights&& weights, int node, int from) {
		for (auto weight : weights) 
			if (weight.node_in1==from && weight.node_out==node) return true;
			else if (weight.node_in2==from && weight.node_out==node) return true;	
		return false;
	}

	template<typename Weights>
	static bool commutes(Weights&& weights, int node1, int node2) {
		return all_of(weights.begin(),weights.end(),
			[node1,node2] (auto& weight) {
				return (node1!=weight.node_in1 || node2!=weight.node_in2) && (node1!=weight.node_in2 || node2!=weight.node_in1);
		});
	}

};




#endif
