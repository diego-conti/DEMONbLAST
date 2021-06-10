#ifndef STORED_COEFFICIENTS_H
#define STORED_COEFFICIENTS_H


#include "partitionprocessor.h"
#include "expressionparser.h"

//stored coefficient lists for a given partition
class StoredCoefficients {
	map<string,CoefficientLists> coefficients_for_diagram;
	bool parse_one(istream& s,const ExpressionParser<StructureConstant>& parser) {
		string diagram;
		getline(s,diagram);
		if (diagram.empty()) return false;
		string line;
		vector<exvector> coefficients;
		while (s && getline(s,line) && !line.empty()) 
			coefficients.push_back(parser.parse_vector_within_brackets(line,","));
		coefficients_for_diagram[diagram]=CoefficientLists(move(coefficients));
		return true;
	}
public:
	StoredCoefficients()=default;
	StoredCoefficients(istream& s) {
		ExpressionParser<StructureConstant> expression_parser;
		while (parse_one(s,expression_parser)) ;
	}
	const CoefficientLists& operator[] (const LabeledTree& diagram) const {
		return coefficients_for_diagram.at(diagram.as_string());
	}
	CoefficientLists& operator[] (const LabeledTree& diagram)  {
		return coefficients_for_diagram[diagram.as_string()];
	}
	void to_stream(ostream& s) const {
		for (auto& pair: coefficients_for_diagram) {
			s<<pair.first<<endl;
			for (auto& v: pair.second)
				s<<v<<endl;
			s<<endl;
		}
	}
};



class PartitionProcessorUsingStoredCoefficients : public PartitionProcessor {
	StoredCoefficients stored_coefficients;
	StoredCoefficients load_stored_coefficients(const vector<int>& partition) {  
	  boost::filesystem::path dir("coefficients");
	  if (!boost::filesystem::is_directory(dir))
  		throw std::runtime_error("directory 'coefficients' does not exist");
		boost::filesystem::path part("coefficients/part"+get_label(partition,"_")+".coeff");
		if (!boost::filesystem::is_regular_file(part))
  		throw std::runtime_error("coefficient file "+part.generic_string()+" not found");
		ifstream s{part.generic_string()};
		return StoredCoefficients{s};
	}	
protected:
	ProcessedDiagram process_single_diagram(LabeledTree& diagram) const override {
		auto& coefficient_lists=stored_coefficients[diagram];
		return processor.process(diagram,coefficient_lists);	
	}
public:
	PartitionProcessorUsingStoredCoefficients(const vector<int>& partition, const DiagramProcessor& processor) 
		: PartitionProcessor(partition,std::move(processor)), stored_coefficients{load_stored_coefficients(partition)} {}
};

class PartitionProcessorStoringCoefficients : public PartitionProcessor {
	mutable StoredCoefficients stored_coefficients;
	void store_coefficients() {  
	  boost::filesystem::path dir("coefficients");
	  if (!boost::filesystem::is_directory(dir)&& !boost::filesystem::create_directories(dir))
	  	throw std::runtime_error("cannot create directory 'coefficients'");	  	
		boost::filesystem::path part("coefficients/part"+get_label(partition,"_")+".coeff");
		ofstream stream{part.generic_string(),std::ofstream::out | std::ofstream::trunc};
  	stored_coefficients.to_stream(stream);
	}	
protected:
	ProcessedDiagram process_single_diagram(LabeledTree& diagram) const override {
		auto& weight_basis=diagram.weight_basis({});
		CoefficientConfigurationWithoutRedundantParameter configuration{WeightBasis{weight_basis}};	
		vector<exvector> coefficients;
		for (auto& group : NiceLieGroup::from_coefficient_configuration(move(configuration)))
			coefficients.push_back(group.c(weight_basis));
		stored_coefficients[diagram]=move(coefficients);
		return PartitionProcessor::process_single_diagram(diagram);
	}
public:
	PartitionProcessorStoringCoefficients(const vector<int>& partition, const DiagramProcessor& processor) 
		: PartitionProcessor(partition,std::move(processor)) {}
	~PartitionProcessorStoringCoefficients() {store_coefficients();}
};


#endif
