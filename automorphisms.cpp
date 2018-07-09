#include "automorphisms.h"

using AssignResult=PartialAutomorphism::AssignResult;

//optimized specialization for labeled trees. TODO optimize also for unlabeled trees.
template<> list<vector<int>> TreeBase<LabeledArrow>::nontrivial_automorphisms() const {
	auto result=automorphisms(*this);
	vector<int> identity(number_of_nodes());
	iota(identity.begin(),identity.end(),0);
	result.remove(identity);
	return result;
}

AssignResult& operator&=(AssignResult& op1, AssignResult op2) {
		if (op1==AssignResult::INCOMPATIBLE || op2==AssignResult::INCOMPATIBLE) 
			op1=AssignResult::INCOMPATIBLE;
		else if (op1==AssignResult::CHANGED || op2==AssignResult::CHANGED) 
			op1=AssignResult::CHANGED; 
		return op1;
}


optional<LabeledArrow> arrow_from_to(int node_in, int node_out, const list<LabeledArrow>& arrows) {
		for (auto& arrow : arrows)
			if (arrow.node_in==node_in && arrow.node_out==node_out) return arrow;
		return nullopt;		
}
optional<LabeledArrow> arrow_from_and(int node_in, int node_in2, const list<LabeledArrow>& arrows) {
		for (auto& arrow : arrows)
			if (arrow.node_in==node_in && arrow.label==node_in2) return arrow;
		return nullopt;		
}	

PartialAutomorphism::AssignResult PartialAutomorphism::match_arrow(LabeledArrow arrow,const list<LabeledArrow>& arrows) {
		AssignResult result=AssignResult::NOCHANGE;
		if (images[arrow.node_in]>=0 && images[arrow.node_out]>=0) {
			auto matching_arrow=arrow_from_to(images[arrow.node_in],images[arrow.node_out],arrows);
			if (matching_arrow) result=assign(arrow.label,matching_arrow.value().label);			
			else return AssignResult::INCOMPATIBLE;
		}
		if (images[arrow.node_in]>=0 && images[arrow.label]>=0) {
			auto matching_arrow=arrow_from_and(images[arrow.node_in],images[arrow.label],arrows);
			if (matching_arrow) result&=assign(arrow.node_out,matching_arrow.value().node_out);
			else return AssignResult::INCOMPATIBLE;
		}		
		return result;
	}

