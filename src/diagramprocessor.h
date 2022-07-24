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

#include "tree.h"
#include "labeled_tree.h"
#include "niceliegroup.h"
#include "niceeinsteinliegroup.h"
#include "partitions.h"
#include "double_arrows_tree.h"
#include "filter.h"
#include "options.h"
#include "diagramanalyzer.h"
#include "coefficientconfiguration.h"

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
  do_not_reorder=64,
  with_lcs_and_ucs=128,
  full_diagram_name=256
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
  virtual ProcessedDiagram process (const LabeledTree& diagram,const CoefficientLists& coefficient_lists) const {
  	//redefined by processor that use the coefficient configuration 
  	throw std::logic_error("fixed coefficient_lists passed to incompatible processor");
  }
  virtual void canonicalize_order(LabeledTree& diagram) const {
  	if (!processing_options().has(ProcessingOption::do_not_reorder))
	  	diagram.canonicalize_order_increasing();
  }
	void set(ProcessingOption option) {processing_options().set(option);}
  void set(DiagramDataOption option) {diagram_data_options().set(option);}
 	void set_sign_configuration_limit(int limit) {diagram_data_options().sign_configurations_limit=limit;}
  void clear(ProcessingOption option) {processing_options().clear(option);}
  void clear(DiagramDataOption option) {diagram_data_options().clear(option);}
  void adapt_to_filter(Filter filter) {diagram_data_options().adapt_to_filter(filter);}
  bool with_automorphisms() const {return processing_options().has(ProcessingOption::with_automorphisms);}
  bool with_diagram_data() const {return diagram_data_options().with_data();}
  bool only_if_lie_algebras() const {return !processing_options().has(ProcessingOption::include_diagrams_no_lie_algebra);}
  bool with_derivations() const {return processing_options().has(ProcessingOption::with_derivations);}
  bool with_polynomial_conditions() const {return processing_options().has(ProcessingOption::with_polynomial_conditions);}
  bool with_enhanced_lie_algebras() const {return processing_options().has(ProcessingOption::with_enhanced_lie_algebras);}
  bool full_diagram_name() const {return processing_options().has(ProcessingOption::full_diagram_name);}
  bool with_lcs_and_ucs() const {return processing_options().has(ProcessingOption::with_lcs_and_ucs);}
  operator DiagramDataOptions() const {return diagram_data_options();}
};

constexpr class with_lie_algebra_tag {} with_lie_algebra;
constexpr class with_nilsoliton_metrics_tag {} with_einstein_metrics;
constexpr class with_ricciflat_metrics_tag {} with_ricciflat_metrics;
constexpr class only_diagrams_tag {} only_diagrams;
constexpr class lie_algebra_table_tag {} lie_algebra_table;

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
  		processor->canonicalize_order(diagram);	//do not reorder if a fixed coefficient_lists was removed
      return processor->process(diagram);
  }
  ProcessedDiagram process(LabeledTree& diagram,const CoefficientLists& coefficient_lists) const {  
      return processor->process(diagram,coefficient_lists);
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
 	void set_sign_configuration_limit(int limit) {processor->set_sign_configuration_limit(limit);}
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
  ProcessedDiagram process(const LabeledTree& diagram,const CoefficientLists& coefficient_lists) const override {
  		return base_processor->process(diagram,coefficient_lists);
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
  ProcessedDiagram process(const LabeledTree& diagram,const CoefficientLists& coefficient_lists) const override {
    ProcessedDiagram processed_by_base = IndirectDiagramProcessor::process(diagram,coefficient_lists);
      DoubleArrowsTree double_tree{diagram};
      processed_by_base.append_extra(double_tree.to_dot_string());
    return processed_by_base;
  }
};

class DiagramProcessorInvertNodes  : public IndirectDiagramProcessor {
	static Weight invert(Weight weight, int nodes) {
		Weight result;
		result.node_in1=nodes-1-weight.node_in2;
		result.node_in2=nodes-1-weight.node_in1;
		result.node_out=nodes-1-weight.node_out;
		return result;
	}
	static bool symmetric(Weight w1, Weight w2, int nodes) {
		int n=nodes-1;
		return w1.node_in1+w2.node_in2 == n && w2.node_in1+w1.node_in2 == n && w1.node_out+w2.node_out == n;
	}
	static exvector invert(int nodes, const exvector& coefficients, const list<Weight>& weights,  const list<Weight>& inverted_weights) {
		exvector result;
		for (auto weight : weights) {
			int i=0;
			for (auto j=inverted_weights.begin();j!=inverted_weights.end();++j,++i) 
				if (symmetric(*j,weight,nodes)) break;
			result.push_back(coefficients[i]);
		}
		return result;
	}
	static CoefficientLists invert(int nodes, const list<Weight>& weights,  const list<Weight>& inverted_weights,const CoefficientLists& coefficient_lists)  {
		vector<exvector> result;
		for (auto& v: coefficient_lists) {
			result.push_back(invert(nodes, v,weights,inverted_weights));
		}
		return CoefficientLists{move(result)};
	}
public:
  using IndirectDiagramProcessor::IndirectDiagramProcessor;
  ProcessedDiagram process(const LabeledTree& diagram) const override {
    LabeledTree inverted_diagram{diagram};    //this discards the WeightBasis object, which has a (small) negative impact on performance
    inverted_diagram.invert_nodes();
    return IndirectDiagramProcessor::process(inverted_diagram); 
  }
	ProcessedDiagram process(const LabeledTree& diagram,const CoefficientLists& coefficient_lists) const override {
    LabeledTree inverted_diagram{diagram};
    inverted_diagram.invert_nodes();
    return IndirectDiagramProcessor::process(inverted_diagram,invert(diagram.number_of_nodes(),diagram.weights(),inverted_diagram.weights(),coefficient_lists)); 
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
  ProcessedDiagram process_diagram_and_configuration(const LabeledTree& diagram,const WeightBasisAndProperties& weight_basis, CoefficientConfiguration&& coefficient_configuration) const {
      auto lie_algebras=NiceLieGroup::from_coefficient_configuration(std::move(coefficient_configuration));
      if (lie_algebras.empty() && only_if_lie_algebras()) return {};
      auto result= process_list(diagram,lie_algebras);
      append_derivations(result,lie_algebras);
 			append_extra(result,diagram);
      if (with_enhanced_lie_algebras()) append_enhanced_lie_algebras(diagram, result,weight_basis);
		  return result;  
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
      CoefficientConfigurationWithoutRedundantParameter configuration{weight_basis};
			return process_diagram_and_configuration(diagram,weight_basis,std::move(configuration));
  }
  ProcessedDiagram process(const LabeledTree& diagram,const CoefficientLists& coefficient_lists) const override {
      auto& weight_basis=diagram.weight_basis(diagram_data_options());
			FixedCoefficientConfiguration configuration{diagram.number_of_nodes(),diagram.weights(),coefficient_lists};
			return process_diagram_and_configuration(diagram,weight_basis,std::move(configuration));
	}
};


vector<int> upper_central_series(const LabeledTree& diagram);
vector<int> lower_central_series(const LabeledTree& diagram);

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
	static string extract_denominators(const exvector& N) {
		stringstream res;
		ex lcm=1;
		ex gcd = N.empty()? 1 : numer(N[0]);
		if (gcd.is_zero()) gcd=1;
		for (auto x: N) {
			lcm=GiNaC::lcm(lcm,denom(x));
			gcd=GiNaC::gcd(gcd,numer(x));
		}
		ex factor=gcd/lcm;
		if (factor!=1) res<<factor;
		res<<"(";
		auto i=N.begin();
		if (i!=N.end()) res<<*i++/factor;
		while (i!=N.end()) res<<","<<*i++/factor;
		res<<")";
		return res.str();
	}

  ProcessedDiagram process_diagram_and_configuration(const LabeledTree& diagram,const WeightBasisAndProperties& weight_basis, CoefficientConfiguration&& coefficient_configuration) const {
      auto groups=NiceLieGroup::from_coefficient_configuration(move(coefficient_configuration));  
      if (groups.empty() && only_if_lie_algebras()) return {};
 	   	string lie_algebras;     
 		  ProgressiveLetter progressive_letter; 
			for (auto group : groups)		{
				string progressive{++progressive_letter};
				if (groups.size()==1) progressive.clear();
			 	lie_algebras+=horizontal(lower_central_series(diagram),"")+":";
			 	lie_algebras+=full_diagram_name()? diagram.name() : diagram.number();
			 	lie_algebras+=progressive+"&"+to_string(group);
			 	if (with_lcs_and_ucs()) 
			 		lie_algebras+="&"+horizontal(upper_central_series(diagram),"");
			  if (with_diagram_data())
				  lie_algebras+="&"+extract_denominators(weight_basis.properties().nikolayevsky_derivation());
				lie_algebras+="&"+weight_basis.properties().classification_of_metrics(group.csquared(weight_basis));
			  lie_algebras+="\\\\\n";
		 }
    	return {lie_algebras,{}};
	}
public:
  ProcessedDiagram process(const LabeledTree& diagram,const CoefficientLists& coefficient_lists) const override {
      auto& weight_basis=diagram.weight_basis(diagram_data_options());
			FixedCoefficientConfiguration configuration{diagram.number_of_nodes(),diagram.weights(),coefficient_lists};
			return process_diagram_and_configuration(diagram,weight_basis,std::move(configuration));
	}
  ProcessedDiagram process(const LabeledTree& diagram) const override {
      auto& weight_basis=diagram.weight_basis(diagram_data_options());
			CoefficientConfigurationWithoutRedundantParameter configuration{WeightBasis{weight_basis}};
			return process_diagram_and_configuration(diagram,weight_basis,std::move(configuration));
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


#endif
