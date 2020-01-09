#ifndef NICEDIAGRAMS_IN_PARTITION_H
#define NICEDIAGRAMS_IN_PARTITION_H

#include "filter.h"
#include "partitions.h"
#include <boost/filesystem.hpp>

class NiceDiagramsInPartition {
	const vector<int> partition;
	list<LabeledTree> trees;	
	NiceDiagramsInPartition(const vector<int>& partition, 	list<LabeledTree>&& trees) : partition{partition},trees{trees} {}
	NiceDiagramsInPartition(const vector<int>& partition, const vector<LabeledTree>& trees) : partition{partition},trees{trees.begin(),trees.end()} {}
public:
	void to_stream(ostream&& s) const {
		to_stream(s);
	}
	void to_stream(ostream& s) const {
		s<<partition.size()<<endl; 
		for (int i : partition) s<<i<<" ";
		s<<endl;
		for (auto& diagram: trees) s<<diagram.as_string()<<endl;	
	}
	static NiceDiagramsInPartition from_stream(istream&& s,const vector<int>& expected_partition) { 
		return from_stream(s,expected_partition);
	}
	static NiceDiagramsInPartition from_stream(istream& s,const vector<int>& expected_partition) {
		int partitionlength;
		s>>partitionlength;
		vector<int> partition(partitionlength);
		for (int i=0;i<partitionlength;++i) s>>partition[i];
		if (expected_partition!=partition) throw std::runtime_error("corrupt diagram data; inconsistent partition");
		list<LabeledTree> trees;	
		int count=0;
		while (s) {
			auto tree=LabeledTree::from_stream(s);
		  if (!tree) break;
		  tree->hash();
	    tree->add_number_to_name(++count);
		  trees.push_back(*tree);
		 }		
		return NiceDiagramsInPartition{partition,move(trees)};
	}
	static NiceDiagramsInPartition compute(const vector<int>& partition, Filter filter={}) {
		int count=0;
		auto diagrams = nice_diagrams(partition,filter,{});
		if (!diagrams.empty()) 
			for (auto& diagram: diagrams)
			  diagram.add_number_to_name(++count);
		return NiceDiagramsInPartition{partition,move(diagrams)};
	}
	
	auto begin() const {return trees.begin();}
	auto end() const {return trees.end();}
	int count() const {return trees.size();}
	void remove_trees(Filter filter, DiagramDataOptions options) {
		::remove_trees(trees, filter,options);
	}
	
	class Iterator {
		list<LabeledTree>::iterator i;
		const DiagramProcessor& processor;
	public:
		Iterator(list<LabeledTree>::iterator i,		const DiagramProcessor& processor) : i{i},processor{processor} {}
		bool operator!=(Iterator it) const {return i!=it.i;}
		ProcessedDiagram process() {return processor.process(*i);}
		auto operator*() {return *i;}
		auto operator++() {++i; return *this;}
	};
	
	auto begin_process(const DiagramProcessor& processor) {
		return Iterator{trees.begin(),processor};	
	}
	auto end_process(const DiagramProcessor& processor) {
		return Iterator{trees.end(),processor};	
	}
};


NiceDiagramsInPartition nice_diagrams_in_partition(const vector<int>& partition) {
  boost::filesystem::path dir("diagrams");
  if (!boost::filesystem::is_directory(dir) && !boost::filesystem::create_directories(dir))
  	throw std::runtime_error("cannot create directory 'diagrams'");
  	
  boost::filesystem::path part("diagrams/part"+get_label(partition,"_")+".diag");
	if (boost::filesystem::is_regular_file(part)) 
		return NiceDiagramsInPartition::from_stream(ifstream{part.generic_string()},partition);
	else {
  	auto result=NiceDiagramsInPartition::compute(partition);
  	result.to_stream(ofstream{part.generic_string(),std::ofstream::out | std::ofstream::trunc});
  	return result;
  }
}
			
NiceDiagramsInPartition nice_diagrams_in_partition(const vector<int>& partition,Filter filter, DiagramDataOptions options) {
	if (filter.has_N1N2N3()) //nonnice diagrams are not cached 
		return NiceDiagramsInPartition::compute(partition,filter);
	auto all_diagrams = nice_diagrams_in_partition(partition);
	all_diagrams.remove_trees(filter,options);
	return all_diagrams;
}	

#endif
