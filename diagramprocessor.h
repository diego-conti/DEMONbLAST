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
#ifndef DIAGRAM_PROCESSOR_H
#define DIAGRAM_PROCESSOR_H

#include <iostream>
#include "tree.h"
#include "labeled_tree.h"
#include "wedge/liesubgroup.h"
#include "niceliegroup.h"
#include "niceeinsteinliegroup.h"
#include "partitions.h"
#include "double_arrows_tree.h"
#include "filter.h"
#include "options.h"
#include "diagramanalyzer.h"

struct ProcessedDiagram {
  string data;
  string extra_data;
  bool empty() {return data.empty();}
  void append_extra(const string& extra) {
    if (extra_data.empty()) extra_data=extra;
    else if (!extra.empty()) extra_data+="\n"+extra;
  }
};

enum class ProcessingOption : unsigned int {
  dflt=0,
  with_automorphisms=2,
  include_diagrams_no_lie_algebra=4,	//this is a filter but it depends on the lie algebras, that are computed at this stage
  with_derivations=8,
  with_polynomial_conditions=16,
  with_enhanced_lie_algebras=32,
  do_not_reorder=64
};



using ProcessingOptions = Options<ProcessingOption>;

class DiagramProcessorImpl {
  ProcessingOptions processing_options_;
  DiagramDataOptions diagram_data_options_;
  friend class IndirectDiagramProcessor;
protected:
  void append_extra(ProcessedDiagram& processed_diagram, const LabeledTree& diagram) const {
			processed_diagram.append_extra(diagram.as_string());
			if (with_diagram_data()) processed_diagram.append_extra(diagram.weight_basis(diagram_data_options()).properties().diagram_data());
      if (with_automorphisms() && !diagram.arrows().empty()) 
        processed_diagram.append_extra(nontrivial_automorphisms_to_string(diagram.nontrivial_automorphisms()));
  }
  template<typename LieAlgebras>
  void append_derivations(ProcessedDiagram& processed_diagram,const LieAlgebras& lie_algebras) const {
  		if (with_derivations())
		  	for (auto& lie_algebra : lie_algebras)
  				 processed_diagram.append_extra(lie_algebra.derivations());
  }
  virtual ProcessingOptions& processing_options() {return processing_options_;} 
  virtual DiagramDataOptions& diagram_data_options() {return diagram_data_options_;}
  virtual ProcessingOptions processing_options() const {return processing_options_;} 
  virtual DiagramDataOptions diagram_data_options() const {return diagram_data_options_;}
public:
  virtual ProcessedDiagram process (const LabeledTree& diagram) const {
			ProcessedDiagram result {diagram.to_dot_string(),{}};
			append_extra(result,diagram);
			return result;
  }
  virtual void canonicalize_order(LabeledTree& diagram) const {
  	if (!processing_options().has(ProcessingOption::do_not_reorder))
	  	diagram.canonicalize_order_increasing();
  }
	void set(ProcessingOption option) {processing_options().set(option);}
  void set(DiagramDataOption option) {diagram_data_options().set(option);}
  void clear(ProcessingOption option) {processing_options().clear(option);}
  void clear(DiagramDataOption option) {diagram_data_options().clear(option);}
  void adapt_to_filter(Filter filter) {diagram_data_options().adapt_to_filter(filter);}
  bool with_automorphisms() const {return processing_options().has(ProcessingOption::with_automorphisms);}
  bool with_diagram_data() const {return diagram_data_options().with_data();}
  bool only_if_lie_algebras() const {return !processing_options().has(ProcessingOption::include_diagrams_no_lie_algebra);}
  bool with_derivations() const {return processing_options().has(ProcessingOption::with_derivations);}
  bool with_polynomial_conditions() const {return processing_options().has(ProcessingOption::with_polynomial_conditions);}
  bool with_enhanced_lie_algebras() const {return processing_options().has(ProcessingOption::with_enhanced_lie_algebras);}
  operator DiagramDataOptions() const {return diagram_data_options();}
};

class with_lie_algebra_tag {} with_lie_algebra;
class with_nilsoliton_metrics_tag {} with_einstein_metrics;
class with_ricciflat_metrics_tag {} with_ricciflat_metrics;
class only_diagrams_tag {} only_diagrams;
class lie_algebra_table_tag {} lie_algebra_table;

class DiagramProcessor  {
  Filter filter_; //REFACTOR: consider removing the filter from this class (impacts nice.cpp)
  unique_ptr<DiagramProcessorImpl> processor;
public:
  DiagramProcessor(only_diagrams_tag) : processor{new DiagramProcessorImpl()} {}
  DiagramProcessor(with_lie_algebra_tag);
  DiagramProcessor(lie_algebra_table_tag);    
  DiagramProcessor(with_nilsoliton_metrics_tag);  
  DiagramProcessor(with_ricciflat_metrics_tag);  
  ProcessedDiagram process(LabeledTree& diagram) const {  
 			processor->canonicalize_order(diagram);
      return processor->process(diagram);
  }
  void invert_nodes();
  void with_delta_otimes_delta();
  void setFilter(Filter newFilter) {
    filter_=newFilter;
    processor->adapt_to_filter(filter_);
  } 
  const Filter& filter() const {return filter_;}
  void set(ProcessingOption option) {processor->set(option);}
  void set(DiagramDataOption option) {processor->set(option);}
  void clear(ProcessingOption option) {processor->clear(option);}
  void clear(DiagramDataOption option) {processor->clear(option);}
  operator DiagramDataOptions() const {return static_cast<DiagramDataOptions>(*processor);}
};

class IndirectDiagramProcessor : public DiagramProcessorImpl {
public:
	IndirectDiagramProcessor(unique_ptr<DiagramProcessorImpl> processor) : base_processor(move(processor)) {}
	ProcessedDiagram process (const LabeledTree& diagram) const override {
		return base_processor->process(diagram);
	}
protected:
  virtual ProcessingOptions& processing_options() override {return base_processor->processing_options();}
  virtual DiagramDataOptions& diagram_data_options() override {return base_processor->diagram_data_options();}
  virtual ProcessingOptions processing_options() const override {return base_processor->processing_options();}
  virtual DiagramDataOptions diagram_data_options() const override {return base_processor->diagram_data_options();}
private:
  unique_ptr<DiagramProcessorImpl> base_processor;
};

class DiagramProcessorWithDeltaOtimesDelta : public IndirectDiagramProcessor {
public:
  using IndirectDiagramProcessor::IndirectDiagramProcessor;
  ProcessedDiagram process(const LabeledTree& diagram) const override {
    ProcessedDiagram processed_by_base = IndirectDiagramProcessor::process(diagram);
      DoubleArrowsTree double_tree{diagram};
      processed_by_base.append_extra(double_tree.to_dot_string());
    return processed_by_base;
  }
};

class DiagramProcessorInvertNodes  : public IndirectDiagramProcessor {
public:
  using IndirectDiagramProcessor::IndirectDiagramProcessor;
  ProcessedDiagram process(const LabeledTree& diagram) const override {
    LabeledTree inverted_diagram{diagram};    //this discards the WeightBasis object, which has a (small) negative impact on performance
    inverted_diagram.invert_nodes();
    return IndirectDiagramProcessor::process(inverted_diagram); 
  }
	void canonicalize_order(LabeledTree& diagram) const override {
  	diagram.canonicalize_order_decreasing();
  }
};

class DiagramProcessorWithLieAlgebras : public DiagramProcessorImpl {
  virtual ProcessedDiagram process_list(const LabeledTree& diagram,const list<NiceLieGroup>& groups) const {  
    string lie_algebras;      
    auto& weight_basis=diagram.weight_basis(diagram_data_options());
		for (auto& group : groups) {
			lie_algebras+=to_string(group)+"\n";
			lie_algebras+=polynomial_equations_for_existence_of_special_metrics(weight_basis, group)+"\n";
			lie_algebras+=weight_basis.properties().classification_of_metrics(group.csquared(weight_basis))+"\n";
		}
    if (groups.empty()) lie_algebras="no Lie algebra";
		return {diagram.to_dot_string(),lie_algebras};
  }
  void append_enhanced_lie_algebras(ProcessedDiagram& processed_diagram, OrderTwoAutomorphism sigma, const EnhancedWeightBasis& weight_basis) const {
     auto lie_algebras = NiceLieGroup::from_weight_basis(weight_basis);
     if (lie_algebras.empty()) return;     
     for (auto& lie_algebra: lie_algebras)
         processed_diagram.append_extra(sigma.to_string()+":"+to_string(lie_algebra));
      processed_diagram.append_extra(horizontal(DiagramAnalyzer{weight_basis.number_of_nodes(),weight_basis.weights_and_coefficients()}.cycles()));
      processed_diagram.append_extra(DiagramAnalyzer{weight_basis.number_of_nodes(),weight_basis.weights_and_coefficients()}.killing());
  }

  void append_enhanced_lie_algebras(const LabeledTree& tree, ProcessedDiagram& processed_diagram, const WeightBasisAndProperties& weight_basis) const {
  	for (auto& sigma : weight_basis.properties().automorphisms_giving_ricci_flat()) 
 			append_enhanced_lie_algebras(processed_diagram,sigma,EnhancedWeightBasis{tree,sigma});
  }
protected:
	string polynomial_equations_for_existence_of_special_metrics(const WeightBasisAndProperties& weight_basis, const NiceLieGroup& group) const {
		if (!with_polynomial_conditions()) return {};
		auto csquared=group.csquared(weight_basis);
		return weight_basis.properties().polynomial_conditions(csquared);
	}
public:
  ProcessedDiagram process(const LabeledTree& diagram) const override {
      auto& weight_basis=diagram.weight_basis(diagram_data_options());
      auto lie_algebras=NiceLieGroup::from_weight_basis(weight_basis);  
      if (lie_algebras.empty() && only_if_lie_algebras()) return {};
      auto result= process_list(diagram,lie_algebras);
      append_derivations(result,lie_algebras);
 			append_extra(result,diagram);
      if (with_enhanced_lie_algebras()) append_enhanced_lie_algebras(diagram, result,weight_basis);
		  return result;
  }
};


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
  ProcessedDiagram process(const LabeledTree& diagram) const override {
      auto& weight_basis=diagram.weight_basis(diagram_data_options());
      auto groups=NiceLieGroup::from_weight_basis(weight_basis);  
      if (groups.empty() && only_if_lie_algebras()) return {};
 	   	string lie_algebras;     
 		  ProgressiveLetter progressive_letter; 
			for (auto group : groups)		{
				string progressive{++progressive_letter};
				if (groups.size()==1) progressive.clear();
			  lie_algebras+=horizontal(lower_central_series(diagram),"")+":"+diagram.name()+progressive+"&"+to_string(group)+"&"+horizontal(upper_central_series(diagram),"");
			  if (with_diagram_data())
				  lie_algebras+="&("+horizontal(weight_basis.properties().nikolayevsky_derivation())+")";
				lie_algebras+="&"+weight_basis.properties().classification_of_metrics(group.csquared(weight_basis));
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
  ProcessedDiagram process(const LabeledTree& diagram) const override {
      auto& weight_basis=diagram.weight_basis(diagram_data_options());
      auto lie_algebras=NiceEinsteinLieGroup::from_weight_basis(weight_basis,metric_type);
      if (lie_algebras.empty() && only_if_lie_algebras()) return {};
      ProcessedDiagram result= process_list(diagram,lie_algebras);
      append_derivations(result,lie_algebras); 			
			append_extra(result,diagram);
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
DiagramProcessor::DiagramProcessor(with_nilsoliton_metrics_tag) : processor{new DiagramProcessorClassifyingMetricLieAlgebras(MetricType::NONFLAT_NILSOLITON)} {} 
DiagramProcessor::DiagramProcessor(with_ricciflat_metrics_tag) : processor{new DiagramProcessorClassifyingMetricLieAlgebras(MetricType::RICCIFLAT)} {} 
DiagramProcessor::DiagramProcessor(lie_algebra_table_tag) : processor{new DiagramProcessorTableOfLieAlgebras()} {} 

#endif
