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
#include "partitions.h"
#include "taskrunner.h"
#include <boost/program_options.hpp>
#include <boost/exception/diagnostic_information.hpp> 
#include "diagramprocessor.h"
#include "nicediagramsinpartition.h"

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
  return "output"+to_string(dimension)+"/part"+get_label(partition,"_")+".dot";
}

string output_path(const vector<int>& partition, const LabeledTree& diagram) {
  int dimension = std::accumulate(partition.begin(),partition.end(),0);
  return "output"+to_string(dimension)+"/graph"+diagram.name()+".dot";
}


int enumerate_nice_diagrams_in_partition(const vector<int>& partition,ostream& stream,const DiagramProcessor& processor) {
    nice_log<<"enumerating diagrams in partition "<<horizontal(partition)<<endl;
			auto diagrams =nice_diagrams_in_partition(partition,processor.filter(),processor);
				for (auto diagram= diagrams.begin_process(processor);diagram!=diagrams.end_process(processor);++diagram) {
			    auto processed=diagram.process();
	        stream<<processed.data;
	        if (!processed.empty() && !processed.extra_data.empty()) stream<<"/*"<<processed.extra_data<<"*/"<<endl;
	        ofstream{output_path(partition,*diagram),std::ofstream::out | std::ofstream::trunc}<<processed.data;
 		    }
		return diagrams.count();
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
  auto diagrams = nice_diagrams(partition,processor.filter(),processor);
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

void process_single_digraph(const po::variables_map& command_line_variables,const DiagramProcessor& processor) {
	auto tree_description=command_line_variables["digraph"].as<string>();
	stringstream s{tree_description};
	auto diagram = LabeledTree::from_stream(s);
	if (!diagram) {cerr<<"no diagram specified"<<endl; return;}
 	auto processed=processor.process(*diagram);
  cout<<processed.data;
	cout<<processed.extra_data;
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
        else if (command_line_variables.count("digraph")) 
					process_single_digraph(command_line_variables,processor);
        else cerr<<"Either --all-partitions, --partition or --digraph must be specified"<<endl;        
}


DiagramProcessor with_options(const po::variables_map& command_line_variables,DiagramProcessor diagram_processor) {
    if (command_line_variables.count("delta-otimes-delta"))
      diagram_processor.with_delta_otimes_delta();
    if (command_line_variables.count("invert"))
      diagram_processor.invert_nodes();
    if (command_line_variables.count("list-diagram-automorphisms")) {
      diagram_processor.set(ProcessingOption::with_automorphisms);
      diagram_processor.set(DiagramDataOption::with_automorphisms);
    }

    if (command_line_variables.count("derivations")) diagram_processor.set(ProcessingOption::with_derivations);
    if (command_line_variables.count("polynomial")) diagram_processor.set(ProcessingOption::with_polynomial_conditions);
    if (command_line_variables.count("enhanced")) diagram_processor.set(ProcessingOption::with_enhanced_lie_algebras);       
    if (command_line_variables.count("matrix-data")) diagram_processor.set(DiagramDataOption::with_matrix_data);
    if (command_line_variables.count("diagonal-ricci-flat-metrics")) diagram_processor.set(DiagramDataOption::with_diagonal_ricci_flat_metrics);
    if (command_line_variables.count("diagonal-nilsoliton-metrics")) diagram_processor.set(DiagramDataOption::with_diagonal_nilsoliton_metrics);
    if (command_line_variables.count("sigma-compatible-ricci-flat-metrics")) diagram_processor.set(DiagramDataOption::with_sigma_compatible_ricci_flat_metrics);       
    if (command_line_variables.count("all-nice-diagrams") || command_line_variables.count("all-diagrams")) diagram_processor.set(ProcessingOption::include_diagrams_no_lie_algebra);
    if (command_line_variables.count("analyze-diagram")) diagram_processor.set(DiagramDataOption::analyze_diagram);
    if (command_line_variables.count("antidiagonal-ricci-flat-sigma")) diagram_processor.set(DiagramDataOption::with_antidiagonal_ricci_flat_sigma);
    if (command_line_variables.count("legacy-weight-order")) diagram_processor.set(ProcessingOption::do_not_reorder);
    
    Filter filter;
    if (command_line_variables.count("only-traceless-derivations")) filter.only_traceless_derivations();
    if (command_line_variables.count("only-MDelta-surjective")) filter.only_MDelta_surjective();
    if (command_line_variables.count("all-diagrams")) filter.N1N2N3();
    if (command_line_variables.count("only-with-nontrivial-automorphisms")) 
      filter.only_nontrivial_automorphisms();
    if (command_line_variables.count("only-with-metric")) 
    	filter.only_with_metric();
    diagram_processor.setFilter(filter);
    return move(diagram_processor);
}

DiagramProcessor create_diagram_processor(const po::variables_map& command_line_variables) {
  auto mode=command_line_variables["mode"].as<string>();
  if (mode=="einstein") {
	  DiagramProcessor diagram_processor{with_einstein_metrics};
	  diagram_processor.set(DiagramDataOption::with_diagonal_nilsoliton_metrics);	//ensure that X is computed
    return diagram_processor;
  }
  if (mode=="ricciflat") {
	  DiagramProcessor diagram_processor{with_ricciflat_metrics};
	  diagram_processor.set(DiagramDataOption::with_diagonal_ricci_flat_metrics); //ensure that X is computed
    return diagram_processor;
 	}
  else if (mode=="diagram")
    return DiagramProcessor{only_diagrams};
  else if (mode=="table")
    return DiagramProcessor{lie_algebra_table};
  else if (mode!="lie-algebra") 
    throw invalid_argument("unrecognized mode");  
  else return DiagramProcessor{with_lie_algebra};
}

void print_help(string program_name, const po::options_description& desc) {
	cout << "Usage: "<< program_name<< " [options]\n";
	cout << desc;
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
            ("digraph", po::value<string>(),
                  "only process diagram indicated by <arg>")
                  
            ("all-nice-diagrams", "in Lie algebra mode, list all nice diagrams [default: only list nice diagrams that correspond to a Lie algebra]") 
            ("all-diagrams", "list all diagrams satisfying N1 N2 N3 [default: only list nice diagrams]; implies all-nice-diagrams")

            ("mode",po::value<string>()->default_value("lie-algebra"), "select operating mode:\n"
                                           "einstein:\t normalize structure constants to obtain einstein metrics with s=1\n"
                                           "ricciflat:\t normalize structure constants to obtain ricci-flat metrics\n"
                                          "diagram:\t only list diagrams\n"
                                          "table:\t list Lie algebras in LaTeX table\n"
                                          "lie-algebra [default]:\t list diagrams and Lie algebras")
            
            ("delta-otimes-delta", "include \\Delta\\otimes\\Delta")                  
            ("list-diagram-automorphisms", "list the automorphisms of each diagram")
            ("invert",  "invert node numbering") 
            ("parallel-mode",  "use multiple threads") 
            ("matrix-data",  "include data depending on the root matrix (rank, etc.)") 
            ("derivations",  "include Lie algebra derivations in output") 
            ("diagonal-ricci-flat-metrics", "include diagonal Ricci-flat metrics")
            ("diagonal-nilsoliton-metrics", "include diagonal nilsoliton metrics")
            ("sigma-compatible-ricci-flat-metrics", "include sigma-compatible Ricci-flat metrics")
            ("polynomial", "include polynomial conditions for the existence of the indicated metrics")
            ("enhanced","include list of sigma-enhanced Lie algebras")         
            ("analyze-diagram","include data depending on diagram combinatorics")
            ("antidiagonal-ricci-flat-sigma","include order two automorphisms inducing antidiagonal ricci-flat metrics")
            ("legacy-weight-order","maintain the weight order coming from the classification algorithm. This option is independent from --invert.")
            
            ("only-traceless-derivations", "exclude diagrams where (1...1) is not in the span of the rows of M_Delta")
            ("only-MDelta-surjective", "only diagrams where M_Delta is surjective")            
            ("only-with-nontrivial-automorphisms", "only diagrams with nontrivial automorphisms)")           
            ("only-with-metric", "only diagrams which potentially admit a metric (conditions H and L)")           
        ;

				try {
        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);
        if (vm.count("help")) print_help(argv[0],desc);
        else if (vm.count("speed-test")) test_speed();
        else 	process(vm,with_options(vm,create_diagram_processor(vm)));
        return 0;
        }
        catch (const boost::program_options::unknown_option& e) {
        		print_help(argv[0],desc);
        		return 1;
        }
}
