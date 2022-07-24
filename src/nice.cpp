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
#include <filesystem>
#include "tree.h"
#include "labeled_tree.h"
#include "wedge/liesubgroup.h"
#include "niceliegroup.h"
#include "partitions.h"
#include "taskrunner.h"
#include <boost/program_options.hpp>
#include <boost/exception/diagnostic_information.hpp> 
#include "diagramprocessor.h"
#include "partitionprocessor.h"


namespace po = boost::program_options;
using namespace GiNaC;
using namespace std;

void create_directory(int n) {
	std::filesystem::create_directory("output");
	
	std::filesystem::path dir("output/"+to_string(n));
	if (std::filesystem::is_directory(dir))
    std::filesystem::remove_all(dir);
  if (!std::filesystem::create_directory(dir))
    cerr<<"cannot create directory"<<endl;
}



int enumerate_nice_diagrams(const list<vector<int>>& partitions,const ProcessorCreator& processor_creator) {
		int count=0;
		for (auto& partition: partitions) {
			auto partition_processor=processor_creator.create(partition);
			count+=partition_processor->process_all();		 	
		}
		return count;
}

void non_parallel_enumerate_nice_diagrams(int dimension,const ProcessorCreator& processor_creator) {
  enumerate_nice_diagrams(partitions(dimension),processor_creator);
}

void parallel_enumerate_nice_diagrams(list<vector<int>> partitions,const ProcessorCreator& processor_creator) {
  if (partitions.empty()) return;
  int dimension=accumulate(partitions.begin()->begin(),partitions.begin()->end(),0);
  auto output_path = [](const vector<int>& partition) {return PartitionProcessor::output_path(partition);};
  auto process_nice_diagrams_in_partition = [&processor_creator](const vector<int>& partition,ostream& stream) {
 		auto partition_processor=processor_creator.create(partition);
		return partition_processor->process_all(stream);
  };
  TaskRunner runner(process_nice_diagrams_in_partition, output_path, partitions);
  runner.run_and_write_to_file();
}


void parallel_enumerate_nice_diagrams(int dimension,const ProcessorCreator& processor_creator) {
  parallel_enumerate_nice_diagrams(partitions(dimension),processor_creator);
}


void test_speed() {
  auto all_partitions = partitions(9);
  list<vector<int>> some_partitions;
  copy_n(all_partitions.begin(),10,back_inserter(some_partitions));
  auto creator=ProcessorCreator::compute_coefficients(DiagramProcessor{with_lie_algebra});
  parallel_enumerate_nice_diagrams(some_partitions,creator);
}


void process_partitions_to_table(int dimension,const ProcessorCreator& processor_creator, bool with_lcs) {
		for (auto& partition: partitions(dimension)) {
			auto partition_processor=processor_creator.create(partition);
		  stringstream output;
		  partition_processor->process_all(output);
    	if (!output.str().empty()) {
    		if (with_lcs) cout<<"&&"<<horizontal(partition,"")<<":\\\\"<<endl;
    		cout<<output.str();	
    	 }
		}
}

void process_all_partitions(const po::variables_map& command_line_variables,const ProcessorCreator& processor_creator) {
  int dimension=command_line_variables["all-partitions"].as<int>();
  if (command_line_variables["mode"].as<string>()=="table") {
     process_partitions_to_table(dimension,processor_creator,command_line_variables.count("lcs-and-ucs"));   
     return;
   }
  create_directory(dimension);
  if (command_line_variables.count("parallel-mode")) 
    parallel_enumerate_nice_diagrams(dimension,processor_creator);
  else 
    non_parallel_enumerate_nice_diagrams(dimension,processor_creator);
}


void process_single_partition(const po::variables_map& command_line_variables,const ProcessorCreator& processor_creator) {
  auto partition= command_line_variables["partition"].as<vector<int>>();
  cout<<"processing partition "<<horizontal(partition)<<endl;
  cout<<enumerate_nice_diagrams({partition},processor_creator)<< " diagrams processed"<<endl;
}

void process_single_digraph(const po::variables_map& command_line_variables,const ProcessorCreator& processor_creator) {
	auto tree_description=command_line_variables["digraph"].as<string>();
	stringstream s{tree_description};
	auto diagram = LabeledTree::from_stream(s);
	if (!diagram) {cerr<<"no diagram specified"<<endl; return;}
	auto partition=lower_central_series(*diagram);
	for (int i=0;i<partition.size()-1;++i) partition[i]-=partition[i+1];
	auto partition_processor=processor_creator.create(partition);
	partition_processor->process(*diagram,cout);
}

void process(const po::variables_map& command_line_variables,const ProcessorCreator& processor_creator) {
        if (command_line_variables.count("all-partitions")) 
          process_all_partitions(command_line_variables,processor_creator);
        else if (command_line_variables.count("partition")) 
          process_single_partition(command_line_variables,processor_creator);
        else if (command_line_variables.count("digraph")) 
					process_single_digraph(command_line_variables,processor_creator);
        else cerr<<"Either --all-partitions, --partition or --digraph must be specified"<<endl;        
}

tribool boolean_value(const po::variables_map& command_line_variables,const string& variable_name) {
	if (!command_line_variables.count(variable_name)) return indeterminate;
	auto value=command_line_variables[variable_name].as<string>();
	if (value=="true" || value=="yes") return true;
	else if (value=="false" ||value=="no") return false;
	throw invalid_argument("unrecognized state: "+value);	
}

ProcessorCreator with_options(const po::variables_map& command_line_variables,DiagramProcessor diagram_processor) {
    if (command_line_variables.count("delta-otimes-delta"))
      diagram_processor.with_delta_otimes_delta();
    if (command_line_variables.count("invert"))
      diagram_processor.invert_nodes();
    if (command_line_variables.count("list-diagram-automorphisms")) {
      diagram_processor.set(ProcessingOption::with_automorphisms);
      diagram_processor.set(DiagramDataOption::with_automorphisms);
      if (command_line_variables.count("do-not-use-automorphisms"))
      	 throw invalid_argument("--list-diagram-automorphisms and --do-not-use-automorphisms are not compatible options");
    }
    

    if (command_line_variables.count("derivations")) diagram_processor.set(ProcessingOption::with_derivations);
    if (command_line_variables.count("polynomial")) diagram_processor.set(ProcessingOption::with_polynomial_conditions);
    if (command_line_variables.count("enhanced")) diagram_processor.set(ProcessingOption::with_enhanced_lie_algebras);       
    if (command_line_variables.count("lcs-and-ucs")) diagram_processor.set(ProcessingOption::with_lcs_and_ucs);       
    if (command_line_variables.count("full-diagram-name")) diagram_processor.set(ProcessingOption::full_diagram_name);       
    if (command_line_variables.count("matrix-data")) diagram_processor.set(DiagramDataOption::with_matrix_data);
    if (command_line_variables.count("diagonal-ricci-flat-metrics")) diagram_processor.set(DiagramDataOption::with_diagonal_ricci_flat_metrics);
    if (command_line_variables.count("diagonal-nilsoliton-metrics")) diagram_processor.set(DiagramDataOption::with_diagonal_nilsoliton_metrics);
    if (command_line_variables.count("sigma-compatible-ricci-flat-metrics")) diagram_processor.set(DiagramDataOption::with_sigma_compatible_ricci_flat_metrics);
    if (command_line_variables.count("only-riemannian-like")) diagram_processor.set(DiagramDataOption::only_riemannian_like_metrics);           
    if (command_line_variables.count("all-nice-diagrams") || command_line_variables.count("all-diagrams")) diagram_processor.set(ProcessingOption::include_diagrams_no_lie_algebra);
    if (command_line_variables.count("analyze-diagram")) diagram_processor.set(DiagramDataOption::analyze_diagram);
    if (command_line_variables.count("antidiagonal-ricci-flat-sigma")) diagram_processor.set(DiagramDataOption::with_antidiagonal_ricci_flat_sigma);
    if (command_line_variables.count("do-not-use-automorphisms")) diagram_processor.set(DiagramDataOption::do_not_use_automorphisms_to_eliminate_signs);
    if (command_line_variables.count("legacy-weight-order")) diagram_processor.set(ProcessingOption::do_not_reorder);
		if (command_line_variables.count("sign-configuration-limit")) diagram_processor.set_sign_configuration_limit(command_line_variables["sign-configuration-limit"].as<int>());
    
    Filter filter;
    if (command_line_variables.count("only-traceless-derivations")) filter.only_traceless_derivations();
    if (command_line_variables.count("kernel-root-matrix-dimension")) filter.only_kernel_MDelta_dimension(command_line_variables["kernel-root-matrix-dimension"].as<string>());
    if (command_line_variables.count("cokernel-root-matrix-dimension")) filter.only_cokernel_MDelta_dimension(command_line_variables["cokernel-root-matrix-dimension"].as<string>());
    if (command_line_variables.count("all-diagrams")) filter.N1N2N3();
    if (command_line_variables.count("only-with-nontrivial-automorphisms")) 
      filter.only_nontrivial_automorphisms();
    if (command_line_variables.count("only-with-metric")) 
    	filter.only_with_metric();
    if (command_line_variables.count("only-with-ad-invariant-metric")) filter.only_passing_obstruction_for_ad_invariant_metric();
    if (command_line_variables.count("irreducible")) filter.only_irreducible();
		filter.simple_nikolayevsky(boolean_value(command_line_variables,"simple-nikolayevsky"));
    diagram_processor.setFilter(filter);

		auto coefficients=command_line_variables["coefficients"].as<string>();
		if (coefficients=="compute") return ProcessorCreator::compute_coefficients(std::move(diagram_processor));
		else if (coefficients=="load") return ProcessorCreator::load_coefficients(std::move(diagram_processor));
		else if (coefficients=="store") return ProcessorCreator::compute_and_store_coefficients(std::move(diagram_processor));
		else return ProcessorCreator::fixed_coefficients(std::move(diagram_processor),coefficients);		
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
//in prospettiva farei che il valore di default è invece usare le tabelle sul disco. 
            ("coefficients",po::value<string>()->default_value("compute"),"select one of the following:\n"
																				"compute:\t The structure constants are computed\n"
																				"store:\t The structure constants are computed and stored to disk\n"
																				"load:\t The structure constants are loaded from disk\n"
																				"a comma-separated vector containing the structure constants (for use with --digraph)")
            ("delta-otimes-delta", "include \\Delta\\otimes\\Delta")                  
            ("list-diagram-automorphisms", "list the automorphisms of each diagram")
            ("lcs-and-ucs", "in table mode, write lower and upper central series dimension")
            ("full-diagram-name", "in table mode, use full diagram name, including the hash code")
            ("do-not-use-automorphisms", "do not compute diagram automorphisms in order to eliminate equivalent families associated to the same diagram; may result in redundant output")
            ("invert",  "invert node numbering") 
            ("parallel-mode",  "use multiple threads") 
            ("matrix-data",  "include data depending on the root matrix (rank, etc.)") 
            ("derivations",  "include Lie algebra derivations in output") 
            ("diagonal-ricci-flat-metrics", "include diagonal Ricci-flat metrics")
            ("diagonal-nilsoliton-metrics", "include diagonal nilsoliton metrics")
            ("sigma-compatible-ricci-flat-metrics", "include sigma-compatible Ricci-flat metrics")
            ("only-riemannian-like", "for diagonal metrics, only include Riemannian metrics and those obtained by changing the signs by an element of the kernel of the mod 2 root matrix")
            ("polynomial", "include polynomial conditions for the existence of the indicated metrics")
            ("enhanced","include list of sigma-enhanced Lie algebras")         
            ("analyze-diagram","include data depending on diagram combinatorics")
            ("antidiagonal-ricci-flat-sigma","include order two automorphisms inducing antidiagonal ricci-flat metrics")
            ("legacy-weight-order","maintain the weight order coming from the classification algorithm. This option is independent from --invert.")
            ("sign-configuration-limit", po::value<int>(), "for each nice diagram, discard sign configurations after the indicated limit")
                        
            ("only-traceless-derivations", "exclude diagrams where (1...1) is not in the span of the rows of M_Delta")
            ("kernel-root-matrix-dimension", po::value<string>(), "filter diagrams where the root matrix has kernel of dimension (=n,<n,>n)")
            ("cokernel-root-matrix-dimension", po::value<string>(), "filter diagrams where the root matrix has cokernel of dimension  (=n,<n,>n)")
            ("only-with-nontrivial-automorphisms", "only diagrams with nontrivial automorphisms)")           
            ("only-with-metric", "only diagrams which potentially admit a metric (conditions H and L)")           
            ("only-with-ad-invariant-metric", "only diagrams which satisfy the necessary condition on generalized lower/upper central series for the existence of an ad-invariant metric")
            ("irreducible", "only connected nice diagrams (i.e. only irreducible nice Lie algebras)")
            ("simple-nikolayevsky", po::value<string>(), "filter diagrams where the Nikolayevsky derivation is simple (i.e. has distinct eigenvalues) or not")
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
