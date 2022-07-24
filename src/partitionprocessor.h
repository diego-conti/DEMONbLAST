/*  Copyright (C) 2018 by Diego Conti, diego.conti@unimib.it      
                                                                     
    This file is part of DEMONbLAST
	                                                                     
    DEMONbLAST is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    DEMONbLAST is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with DEMONbLAST.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef PARTITION_PROCESSOR_H
#define PARTITION_PROCESSOR_H

#include "nicediagramsinpartition.h"

class PartitionProcessor {
protected:
	const vector<int> partition;
	const DiagramProcessor& processor;
	virtual ProcessedDiagram process_single_diagram(LabeledTree& diagram) const {
		return processor.process(diagram);	
	}
public:
	PartitionProcessor(const vector<int>& partition, const DiagramProcessor& processor) : partition{partition},processor{processor} {}
	void process(LabeledTree& diagram,ostream& s) const {
		auto processed=process_single_diagram(diagram);
		s<<processed.data;
		s<<processed.extra_data<<endl;	
	}	
	int process_all(ostream& s) const {
		auto diagrams =nice_diagrams_in_partition(partition,processor.filter(),processor);
		for (auto diagram : diagrams) {
	   auto processed=process_single_diagram(diagram);
     s<<processed.data;     
     if (!processed.empty() && !processed.extra_data.empty()) s<<"/*"<<processed.extra_data<<"*/"<<endl;
     ofstream{output_path(partition,diagram),std::ofstream::out | std::ofstream::trunc}<<processed.data;
		}
		return diagrams.count();
	}	
	int process_all() const {
		stringstream output;
		auto count=process_all(output);
		if (!output.str().empty())
 			ofstream{output_path(partition),std::ofstream::out | std::ofstream::trunc}<<output.str();
 		return count;
  }
	static string output_path(const vector<int>& partition, const string& filename) {
  	int dimension = std::accumulate(partition.begin(),partition.end(),0);
	  return "output/"+std::to_string(dimension)+"/"+filename;
	}
	static string output_path(const vector<int>& partition) {
  	int dimension = std::accumulate(partition.begin(),partition.end(),0);
	  return "output/"+std::to_string(dimension)+"/part"+get_label(partition,"_")+".dot";
	}
	static string output_path(const vector<int>& partition, const LabeledTree& diagram) {
	  int dimension = std::accumulate(partition.begin(),partition.end(),0);
  	return "output/"+std::to_string(dimension)+"/graph"+diagram.name()+".dot";
	}
	virtual ~PartitionProcessor()=default;
};

enum class ProcessorCreatingMode {
	COMPUTE, COMPUTE_STORE, LOAD, FIXED
};

class ProcessorCreator {
	DiagramProcessor processor;
	ProcessorCreatingMode mode;
	string coefficients;
	ProcessorCreator(DiagramProcessor&& processor, ProcessorCreatingMode mode) : processor{std::move(processor)}, mode{mode} {}
	ProcessorCreator(DiagramProcessor&& processor, ProcessorCreatingMode mode, const string& coefficients) : processor{std::move(processor)}, mode{mode}, coefficients{coefficients} {} 
public:
	static ProcessorCreator compute_coefficients(DiagramProcessor&& processor){return ProcessorCreator(std::move(processor),ProcessorCreatingMode::COMPUTE);}
	static ProcessorCreator compute_and_store_coefficients(DiagramProcessor&& processor){return ProcessorCreator(std::move(processor),ProcessorCreatingMode::COMPUTE_STORE);}
	static ProcessorCreator load_coefficients(DiagramProcessor&& processor) {return ProcessorCreator(std::move(processor),ProcessorCreatingMode::LOAD);}
	static ProcessorCreator fixed_coefficients(DiagramProcessor&& processor, const string& coefficients) {return ProcessorCreator(std::move(processor),ProcessorCreatingMode::FIXED,coefficients);}
	unique_ptr<PartitionProcessor> create(const vector<int>& partition) const;
};

#endif
