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
#include "tree.h"
#include "labeled_tree.h"
#include "wedge/liesubgroup.h"
#include "niceliegroup.h"
#include "niceeinsteinliegroup.h"
#include "partitions.h"
#include "double_arrows_tree.h"
#include "filter.h"

struct ProcessedDiagram {
  string data;
  string extra_data;
  bool empty() {return data.empty();}
  void append_extra(const string& extra) {
    if (extra_data.empty()) extra_data=extra;
    else if (!extra.empty()) extra_data+="\n"+extra;
  }
};

enum class Option : unsigned int {
  dflt=0,
  with_diagram_data=1,
  with_automorphisms=4,
  include_diagrams_no_lie_algebra=8,
  with_derivations=16,
  with_ricci_flat_metrics=32
};

auto operator^=(Option& options, Option option) {
  using underlying=std::underlying_type_t<Option>;
  underlying value=static_cast<underlying>(options)^static_cast<underlying>(option);
  return options= static_cast<Option>(value);
}
auto operator|=(Option& options, Option option) {
  using underlying=std::underlying_type_t<Option>;
  underlying value=static_cast<underlying>(options)|static_cast<underlying>(option);
  return options= static_cast<Option>(value);
}
bool operator&(Option options, Option option) {
  using underlying=std::underlying_type_t<Option>;
  return static_cast<underlying>(options)&static_cast<underlying>(option);
}
            
struct Options  {
  Option options=Option::dflt;
public:
  bool with_automorphisms() const {return options & Option::with_automorphisms;}
  bool with_diagram_data() const {return options & Option::with_diagram_data;}
  bool only_if_lie_algebras() const {return !(options & Option::include_diagrams_no_lie_algebra);}
  bool with_derivations() const {return options & Option::with_derivations;}
  bool with_diagonal_ricci_flat_metrics() const {return options & Option::with_ricci_flat_metrics;}
  void log() const {nice_log<<"options = "<<static_cast<unsigned int>(options)<<endl;}
  void set(Option option) {options|=option;}
  void clear(Option option) {set(option); options^=option;}
};


class DiagramProcessorImpl {
protected:
  void append_extra(ProcessedDiagram& processed_diagram, const LabeledTree& diagram, Options options) const {
			if (options.with_diagram_data()) processed_diagram.append_extra(diagram.weight_basis().properties().diagram_data());  
      if (options.with_automorphisms() && !processed_diagram.empty() && !diagram.arrows().empty()) 
        processed_diagram.append_extra(nontrivial_automorphisms_to_string(diagram.nontrivial_automorphisms()));
  }
  template<typename LieAlgebras>
  void append_derivations(ProcessedDiagram& processed_diagram,const LieAlgebras& lie_algebras,Options options) const {
  		if (options.with_derivations())
		  	for (auto& lie_algebra : lie_algebras)
  				 processed_diagram.append_extra(lie_algebra.derivations());
  }
public:
  virtual ProcessedDiagram process (const LabeledTree& diagram, Options options) const {
			ProcessedDiagram result {diagram.to_dot_string(),{}};
			append_extra(result,diagram,options);
			return result;
  }
};


class with_lie_algebra_tag {} with_lie_algebra;
class with_einstein_metrics_tag {} with_einstein_metrics;
class with_ricciflat_metrics_tag {} with_ricciflat_metrics;
class only_diagrams_tag {} only_diagrams;
class lie_algebra_table_tag {} lie_algebra_table;


class DiagramProcessor  {
  Filter filter_; //REFACTOR: consider removing the filter from this class (impacts nice.cpp)
  Options options;
  unique_ptr<DiagramProcessorImpl> processor;
public:
  DiagramProcessor(only_diagrams_tag) : processor{new DiagramProcessorImpl()} {}
  DiagramProcessor(with_lie_algebra_tag);
  DiagramProcessor(lie_algebra_table_tag);    
  DiagramProcessor(with_einstein_metrics_tag);  
  DiagramProcessor(with_ricciflat_metrics_tag);  
  ProcessedDiagram process(const LabeledTree& diagram) const {  
      options.log();
      return processor->process(diagram,options);
  }
  void invert_nodes();
  void with_delta_otimes_delta();
  void setFilter(Filter newFilter) {
    filter_=newFilter;
  } 
  const Filter& filter() const {return filter_;}
  void set(Option option) {options.set(option);}
  void clear(Option option) {options.clear(option);}
};


class DiagramProcessorWithDeltaOtimesDelta : public DiagramProcessorImpl {
  unique_ptr<DiagramProcessorImpl> base_processor;
public:
  DiagramProcessorWithDeltaOtimesDelta(  unique_ptr<DiagramProcessorImpl> processor) : base_processor(move(processor)) {}
  ProcessedDiagram process(const LabeledTree& diagram, Options options) const override {
    ProcessedDiagram processed_by_base = base_processor->process(diagram,options);
    if (!processed_by_base.empty()) {
      DoubleArrowsTree double_tree{diagram};
      processed_by_base.append_extra(double_tree.to_dot_string());
    }
    return processed_by_base;
  }
};

class DiagramProcessorInvertNodes  : public DiagramProcessorImpl {
  unique_ptr<DiagramProcessorImpl> base_processor;
public:
  DiagramProcessorInvertNodes(  unique_ptr<DiagramProcessorImpl> processor) : base_processor(move(processor)) {}
  ProcessedDiagram process(const LabeledTree& diagram, Options options) const override {
    LabeledTree inverted_diagram{diagram};    //this discards the WeightBasis object, which has a (small) negative impact on performance
    inverted_diagram.invert_nodes();
    return base_processor->process(inverted_diagram,options); 
  }
};



class DiagramProcessorWithLieAlgebras : public DiagramProcessorImpl {
  virtual ProcessedDiagram process_list(const LabeledTree& diagram,const list<NiceLieGroup>& groups, Options options) const {  
    string lie_algebras;      
		for (auto& group : groups) {
			lie_algebras+=to_string(group)+"\n";
			if (options.with_diagonal_ricci_flat_metrics()) lie_algebras+=polynomial_equations_for_existence_of_ricci_flat_metric(diagram.weight_basis(), group);
		}
    if (groups.empty()) lie_algebras="no Lie algebra";
		return {diagram.to_dot_string(),lie_algebras};
  }
  ProcessedDiagram process_list_and_automorphisms(const LabeledTree& diagram,const list<NiceLieGroup>& groups, Options options) const {  
    auto result=process_list(diagram,groups,options);
    if (groups.size()<2) return result; //no need to apply automorphisms, since we only have one Lie algebra
    auto automorphisms= diagram.nontrivial_automorphisms();
    if (automorphisms.empty()) return result;
    result.append_extra("More than one Lie algebra and more than one automorphism:");
 		for (auto automorphism : automorphisms) {
 		  result.append_extra("sigma = "+permutation_to_string(automorphism));
  		for (auto group : groups) 
	  	  result.append_extra(to_string(change_basis(group,automorphism)));
    }
    return result;
  }
protected:
	static string polynomial_equations_for_existence_of_ricci_flat_metric(const WeightBasis& weight_basis, const NiceLieGroup& group) {
		auto csquared=group.csquared(weight_basis);
		auto equations=weight_basis.properties().polynomial_equations_for_existence_of_ricci_flat_metric(csquared);
		return equations.empty()? string{} : "Ricci-flat when "+horizontal(equations)+", X="+horizontal(weight_basis.properties().ricci_flat_X_ijk());
	}
public:
  ProcessedDiagram process(const LabeledTree& diagram, Options options) const override {
      auto& weight_basis=diagram.weight_basis();
      auto lie_algebras=NiceLieGroup::from_weight_basis(weight_basis);  
      if (lie_algebras.empty() && options.only_if_lie_algebras()) return {};        
      auto result= options.with_automorphisms() ? process_list_and_automorphisms(diagram,lie_algebras,options) : process_list(diagram,lie_algebras,options);
      append_derivations(result,lie_algebras,options);
 			append_extra(result,diagram,options);
		  return result;
  }
};


set<int> consecutive_numbers(int begin, int end) {
  set<int> numbers;
  while (begin!=end) numbers.emplace(begin++);
  return numbers;
}

vector<int> upper_central_series(const LabeledTree& diagram,set<int> subspace) {
  set<int> nodes=consecutive_numbers(0,diagram.number_of_nodes());
  for (auto& arrow: diagram.arrows()) 
      if (subspace.find(arrow.node_out)==subspace.end()) {nodes.erase(arrow.node_in);}
  for (int node: nodes) subspace.emplace(node);
  if (subspace.size()==diagram.number_of_nodes()) return {diagram.number_of_nodes()};
  auto rest_of_series= upper_central_series(diagram, subspace);
  rest_of_series.insert(rest_of_series.begin(),nodes.size());
  return rest_of_series;
}

vector<int> upper_central_series(const LabeledTree& diagram) {
  return upper_central_series(diagram, {});
}

vector<int> lower_central_series(const LabeledTree& diagram, set<int> subspace) {
  if (subspace.empty()) return {};
  set<int> next_subspace;
  for (auto& arrow: diagram.arrows()) 
      if (subspace.find(arrow.node_in)!=subspace.end()) next_subspace.insert(arrow.node_out);
  auto rest_of_series= lower_central_series(diagram, next_subspace);
  rest_of_series.insert(rest_of_series.begin(),subspace.size());
  return rest_of_series;
}

vector<int> lower_central_series(const LabeledTree& diagram) {
  return lower_central_series(diagram, consecutive_numbers(0,diagram.number_of_nodes()));
}


class ProgressiveLetter {
	string representation;
public:
	ProgressiveLetter& operator++() {
		int increase_at=representation.size()-1;
		while (increase_at>=0 && ++representation.at(increase_at)>'z') {
			representation.at(increase_at)='a';
			--increase_at;
		}
		if (increase_at<0) representation.push_back('a');
		return *this;
	}
	operator string() const {return representation;}
};
class DiagramProcessorTableOfLieAlgebras : public DiagramProcessorWithLieAlgebras {
public:
  ProcessedDiagram process_list(const LabeledTree& diagram,const list<NiceLieGroup>& groups, Options options) const override {  
    string lie_algebras;     
    ProgressiveLetter progressive_letter; 
		for (auto group : groups)		{
			string progressive{++progressive_letter};
			if (groups.size()==1) progressive.clear();
		  lie_algebras+=horizontal(lower_central_series(diagram),"")+":"+diagram.name()+progressive+"&"+to_string(group)+"&"+horizontal(upper_central_series(diagram),"");
 			if (options.with_diagonal_ricci_flat_metrics()) lie_algebras+="&"+polynomial_equations_for_existence_of_ricci_flat_metric(diagram.weight_basis(), group);
		  lie_algebras+="\\\\\n";
		 }
    return {lie_algebras,{}};
  }
};

class DiagramProcessorClassifyingMetricLieAlgebras : public DiagramProcessorImpl {
  ProcessedDiagram process_list(const LabeledTree& diagram,const list<NiceEinsteinLieGroup>& groups) const {
      string lie_algebras;      
			for (auto group : groups) 
				lie_algebras+=to_string(group)+"\n";
        if (groups.empty()) lie_algebras="no Lie algebra";	  
      return {diagram.to_dot_string(),lie_algebras};
  }
  MetricType metric_type;
public:
	DiagramProcessorClassifyingMetricLieAlgebras(MetricType metric_type) : metric_type{metric_type} {}
  ProcessedDiagram process(const LabeledTree& diagram, Options options) const override {
      auto& weight_basis=diagram.weight_basis();
      auto lie_algebras=NiceEinsteinLieGroup::from_weight_basis(weight_basis,metric_type);          
      if (lie_algebras.empty() && options.only_if_lie_algebras()) return {};        
      ProcessedDiagram result= process_list(diagram,lie_algebras);
      append_derivations(result,lie_algebras,options); 			
			append_extra(result,diagram,options);
			return result;
		}
};



void DiagramProcessor::invert_nodes() {
    if (dynamic_cast<DiagramProcessorInvertNodes*>(processor.get())==nullptr)
    processor =  make_unique<DiagramProcessorInvertNodes>(move(processor));
}

void DiagramProcessor::with_delta_otimes_delta() {
    if (dynamic_cast<DiagramProcessorWithDeltaOtimesDelta*>(processor.get())==nullptr)
     processor =  make_unique<DiagramProcessorWithDeltaOtimesDelta>(move(processor));
}


DiagramProcessor::DiagramProcessor(with_lie_algebra_tag) : processor{new DiagramProcessorWithLieAlgebras()} {}
DiagramProcessor::DiagramProcessor(with_einstein_metrics_tag) : processor{new DiagramProcessorClassifyingMetricLieAlgebras(MetricType::NONFLAT_EINSTEIN)} {} 
DiagramProcessor::DiagramProcessor(with_ricciflat_metrics_tag) : processor{new DiagramProcessorClassifyingMetricLieAlgebras(MetricType::RICCIFLAT)} {} 
DiagramProcessor::DiagramProcessor(lie_algebra_table_tag) : processor{new DiagramProcessorTableOfLieAlgebras()} {} 

