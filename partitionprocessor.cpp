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

#include "partitionprocessor.h"
#include "storedcoefficients.h"
#include "expressionparser.h"

class PartitionProcessorUsingFixedCoefficients : public PartitionProcessor {
	CoefficientLists coefficients;
	
	static CoefficientLists create_coefficients(const string& coefficients) {
		ExpressionParser<StructureConstant> expression_parser;
		exvector v=expression_parser.parse_vector(coefficients,",");
		return {v};
	}
	
protected:
	ProcessedDiagram process_single_diagram(LabeledTree& diagram) const override {
		return processor.process(diagram,coefficients);	
	}
public:
	PartitionProcessorUsingFixedCoefficients(const vector<int>& partition, const DiagramProcessor& processor, const string& coefficients)
		: PartitionProcessor(partition,processor), coefficients{create_coefficients(coefficients)}  {}
};


unique_ptr<PartitionProcessor> ProcessorCreator::create(const vector<int>& partition) const {
	switch (mode) {
		case ProcessorCreatingMode::COMPUTE :
			return make_unique<PartitionProcessor>(partition,processor);
		case ProcessorCreatingMode::COMPUTE_STORE :
			return make_unique<PartitionProcessorStoringCoefficients>(partition,processor);
		case ProcessorCreatingMode::LOAD :
			return make_unique<PartitionProcessorUsingStoredCoefficients>(partition,processor);
		case ProcessorCreatingMode::FIXED :			
			return make_unique<PartitionProcessorUsingFixedCoefficients>(partition,processor,coefficients);
		default:
			return nullptr;
		}
}

