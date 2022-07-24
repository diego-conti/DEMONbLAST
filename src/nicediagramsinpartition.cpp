#include "nicediagramsinpartition.h"

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
