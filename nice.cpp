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
#include <iostream>
#include <fstream>
#include <ginac/ginac.h>
#include <future>
#include <boost/filesystem.hpp>
#include "tree.h"
#include "labeled_tree.h"
#include "wedge/liesubgroup.h"
#include "niceliegroup.h"
#include "niceeinsteinliegroup.h"
#include "partitions.h"
#include "taskrunner.h"
#include <boost/program_options.hpp>
#include "diagramprocessor.h"


namespace po = boost::program_options;
using namespace GiNaC;

void create_directory(int n) {
  boost::filesystem::path dir("output"+to_string(n));
  if (boost::filesystem::is_directory(dir)) 
    boost::filesystem::remove_all(dir);
  if (!boost::filesystem::create_directories(dir))
    cerr<<"cannot create directory"<<endl;
}

string output_path(const vector<int>& partition, string filename) {
  int dimension = std::accumulate(partition.begin(),partition.end(),0);
  return "output"+to_string(dimension)+"/"+filename;
}
string output_path(const vector<int>& partition) {
  int dimension = std::accumulate(partition.begin(),partition.end(),0);
  return "output"+to_string(dimension)+"/graph"+get_label(partition,"_")+".dot";
}



int enumerate_nice_diagrams_in_partition(const vector<int>& partition,ostream& stream,const DiagramProcessor& processor) {
    nice_log<<"enumerating diagrams in partition "<<horizontal(partition)<<endl;
      map<int,int> hashes_to_hashes;
      int count=0;
			auto diagrams = nice_diagrams(partition,processor.filter());
			if (!diagrams.empty()) {
				for (auto& diagram: diagrams) {
				  ++hashes_to_hashes[diagram.hash()];
				  diagram.add_number_to_name(++count);
				  nice_log<<diagram.name()<<endl;
			    auto processed=processor.process(diagram);
	        stream<<processed.data;
	        if (!processed.empty() && !processed.extra_data.empty()) stream<<"/*"<<processed.extra_data<<"*/"<<endl;
 		    }
  		}
		for (auto p : hashes_to_hashes) 
		  if (p.second>1) 
		    nice_log<<" multiple ("<<p.second<<") diagrams with hash "<<p.first<<endl;		  
		return count;
}

int enumerate_nice_diagrams(const list<vector<int>>& partitions,const DiagramProcessor& processor) {
		int count=0;
		for (auto& partition: partitions) {
		  stringstream output;
			count+=enumerate_nice_diagrams_in_partition(partition,output,processor);
			if (!output.str().empty())
  			ofstream{output_path(partition),std::ofstream::out | std::ofstream::trunc}<<output.str();
		}
		return count;
}

void process_single_diagram(const vector<int>& partition, int hash,const DiagramProcessor& processor) {
  int dimension = std::accumulate(partition.begin(),partition.end(),0);
  auto diagrams = nice_diagrams(partition,processor.filter());
  bool found=false;
	for (auto& diagram : diagrams) 
	  if (diagram.hash()==hash) {
	        found=true;
	        auto processed=processor.process(diagram);
	        cout<<processed.data;
	        cout<<processed.extra_data;
		}
  if (!found)
  	cout<<"diagram not found"<<endl;  
}

void non_parallel_enumerate_nice_diagrams(int dimension,const DiagramProcessor& processor) {
  enumerate_nice_diagrams(partitions(dimension),processor);
}

void parallel_enumerate_nice_diagrams(list<vector<int>> partitions,const DiagramProcessor& processor) {
  if (partitions.empty()) return;
  int dimension=accumulate(partitions.begin()->begin(),partitions.begin()->end(),0);
  string (&output_path) (const vector<int>&) = ::output_path;
  auto process_nice_diagrams_in_partition = [&processor](const vector<int>& partition,ostream& stream) {
     return enumerate_nice_diagrams_in_partition(partition,stream,processor);
  };
  TaskRunner runner(process_nice_diagrams_in_partition, output_path, partitions);
  runner.run_and_write_to_file();
}


void parallel_enumerate_nice_diagrams(int dimension,const DiagramProcessor& processor) {
  parallel_enumerate_nice_diagrams(partitions(dimension),processor);
}




void test_speed() {
  auto all_partitions = partitions(9);
  list<vector<int>> some_partitions;
  copy_n(all_partitions.begin(),10,back_inserter(some_partitions));
  parallel_enumerate_nice_diagrams(some_partitions,DiagramProcessor{with_lie_algebra});
}


void process_partitions_to_table(int dimension,const DiagramProcessor& processor) {
		for (auto& partition: partitions(dimension)) { 
		  stringstream output;
    	enumerate_nice_diagrams_in_partition(partition,output,processor);	
    	if (!output.str().empty()) {
     	 cout<<"&&"<<horizontal(partition,"")<<":\\\\"<<endl<<output.str();
    	 }
		}
}

void process_all_partitions(const po::variables_map& command_line_variables,const DiagramProcessor& processor) {
  int dimension=command_line_variables["all-partitions"].as<int>();
  if (command_line_variables["mode"].as<string>()=="table") {
     process_partitions_to_table(dimension,processor);   
     return;
   }
  create_directory(dimension);
  if (command_line_variables.count("parallel-mode")) 
    parallel_enumerate_nice_diagrams(dimension,processor);
  else 
    non_parallel_enumerate_nice_diagrams(dimension,processor);
}


void process_single_partition(const po::variables_map& command_line_variables,const DiagramProcessor& processor) {
  auto partition= command_line_variables["partition"].as<vector<int>>();
  if (command_line_variables.count("diagram"))  {
       int diagram_hash=command_line_variables["diagram"].as<int>();
       cout<<"processing diagram(s) with hash "<<diagram_hash<<endl<<" partition "<<horizontal(partition)<<endl;
       process_single_diagram(partition,diagram_hash,processor);
   }
   else if (command_line_variables.count("partition"))
    {
      cout<<"processing partition "<<horizontal(partition)<<endl;
      ofstream file{output_path(partition),std::ofstream::out | std::ofstream::trunc};
      cout<<enumerate_nice_diagrams_in_partition(partition,file,processor)<< " diagrams processed"<<endl;
    }    
}

void process(const po::variables_map& command_line_variables,const DiagramProcessor& processor) {
        if (command_line_variables.count("all-partitions")) 
          process_all_partitions(command_line_variables,processor);
        else if (command_line_variables.count("partition")) 
          process_single_partition(command_line_variables,processor);
        else cerr<<"Either --all-partitions or --partition must be specified"<<endl;        
}


DiagramProcessor with_options(const po::variables_map& command_line_variables,DiagramProcessor diagram_processor) {
    if (command_line_variables.count("delta-otimes-delta"))
      diagram_processor.with_delta_otimes_delta();
    if (command_line_variables.count("invert"))
      diagram_processor.invert_nodes();
    if (command_line_variables.count("list-diagram-automorphisms")) diagram_processor.set(Option::with_automorphisms);
    if (command_line_variables.count("diagram-data")) diagram_processor.set(Option::with_diagram_data);
    if (command_line_variables.count("derivations")) diagram_processor.set(Option::with_derivations);
    if (command_line_variables.count("all-nice-diagrams") || command_line_variables.count("all-diagrams")) diagram_processor.set(Option::include_diagrams_no_lie_algebra);
    Filter filter;
    if (command_line_variables.count("no-coordinate-hyperplane")) filter.no_coordinate_hyperplane();
    if (command_line_variables.count("only-traceless-derivations")) filter.only_traceless_derivations();
    if (command_line_variables.count("only-MDelta-surjective")) filter.only_MDelta_surjective();
    if (command_line_variables.count("only-reachable-orthant")) filter.only_reachable_orthant();
    if (command_line_variables.count("all-diagrams")) filter.N1N2N3();
    diagram_processor.setFilter(filter);     
    return move(diagram_processor);
}

DiagramProcessor create_diagram_processor(const po::variables_map& command_line_variables) {
  auto mode=command_line_variables["mode"].as<string>();
  if (mode=="einstein")
    return DiagramProcessor{with_einstein_metrics};
  else if (mode=="diagram")
    return DiagramProcessor{only_diagrams};
  else if (mode=="table")
    return DiagramProcessor{lie_algebra_table};
  else if (mode!="lie-algebra") 
    throw invalid_argument("unrecognized mode");  
  else return DiagramProcessor{with_lie_algebra};
}


int main(int argc, char* argv[]) {
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help", "produce help message")
            ("speed-test", "run a speed test")
            
            ("all-partitions", po::value<int>(),
                  "all partitions of dimension <arg>")
            ("partition", po::value<vector<int>>()->multitoken(),
                  "only process partition arg")
            ("diagram", po::value<int>(),
                  "only process diagram with hash <arg>")
                  
            ("all-nice-diagrams", "in Lie algebra mode, list all nice diagrams [default: only list nice diagrams that correspond to a Lie algebra]") 
            ("all-diagrams", "list all diagrams satisfying N1 N2 N3 [default: only list nice diagrams]; implies all-nice-diagrams")

            ("mode",po::value<string>()->default_value("lie-algebra"), "select operating mode:\n"
                                           "einstein:\t discard diagrams with M_Delta not surjective; normalize structure constants to obtain einstein metric\n"
                                          "diagram:\t only list diagrams\n"
                                          "table:\t list Lie algebras in LaTeX table\n"
                                          "lie-algebra [default]:\t list diagrams and Lie algebras")
            
            ("delta-otimes-delta", "include \\Delta\\otimes\\Delta")                  
            ("list-diagram-automorphisms", "list the automorphisms of each diagram")
            ("invert",  "invert node numbering") 
            ("parallel-mode",  "use multiple threads") 
            ("diagram-data",  "include diagram data in output") 
            ("derivations",  "include Lie algebra derivations in output") 

            
            ("only-traceless-derivations", "exclude diagrams where (1...1) is not in the span of the rows of M_Delta")
            ("no-coordinate-hyperplane", "exclude diagrams where X_ijk is in coordinate hyperplane (implies --only-traceless-derivations)")
            ("only-MDelta-surjective", "only diagrams where M_Delta is surjective")            
            ("only-reachable-orthant", "only diagrams where X_ijk intersects reachable orthant (implies --only-traceless-derivations; use with --only-MDelta-surjective)")
            
        ;

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);
        if (vm.count("help")) {
            cout << "Usage: "<< argv[0]<< " [options]\n";
            cout << desc;
        }        
        else if (vm.count("speed-test")) {
        	test_speed();
        	return 0;
        }
        else process(vm,with_options(vm,create_diagram_processor(vm)));
        return 0;
}
