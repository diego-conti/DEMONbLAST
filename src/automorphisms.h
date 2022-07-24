#ifndef AUTOMORPHISMS_H
#define AUTOMORPHISMS_H
#include "tree.h"


class PartialAutomorphism {
public:
	enum class AssignResult {NOCHANGE, CHANGED, INCOMPATIBLE};	

	PartialAutomorphism(const vector<int>& hash_codes ) : images(hash_codes.size(),UNASSIGNED), hash_codes_of_unassigned{hash_codes} {
		for (auto& code : hash_codes_of_unassigned)
			if (code==IN_IMAGE) code=-2;	//IN_IMAGE value is reserved, so replace occurrences of IN_IMAGE with a different value.
	}
	
	
	template<typename ArrowType>
	list<vector<int>> enlargements(const TreeBase<ArrowType>& tree) const {
		if (complete_permutation()) 
			return {images};
		list<vector<int>> result;
		int i=first_unassigned();
		for (int j=0;j<images.size();++j) 
			if (hash_codes_of_unassigned[j]==tree.node_hash()[i]) {
				auto enlarged=with_assignment(tree.arrows(), i,j);
				if (enlarged) result.splice(result.end(),enlarged.value().enlargements(tree));
			}
		return result;
	}
private:

	static constexpr int UNASSIGNED=-1;
	static constexpr int IN_IMAGE=-1;
	vector<int> images;											//i-th element is either sigma_i or UNASSIGNED
	vector<int> hash_codes_of_unassigned;		//i-th  element is either the hash chode of i-th node or IN_IMAGE if i=image[j] for some j
	
	AssignResult assign(int i, int sigma_i) {
		if (images[i]==sigma_i) return AssignResult::NOCHANGE;
		else if (images[i]==UNASSIGNED && hash_codes_of_unassigned[sigma_i]!=IN_IMAGE) {
			images[i]=sigma_i;
			hash_codes_of_unassigned[sigma_i]=IN_IMAGE;
			return AssignResult::CHANGED;
		}
		else return AssignResult::INCOMPATIBLE;
	}
	
//TODO for unlabeled arrows, the logic should be that assignments are made when only one choice is possible
	AssignResult match_arrow(LabeledArrow arrow,const list<LabeledArrow>& arrows);

	template<typename ArrowType>
	optional<PartialAutomorphism> with_assignment(const list<ArrowType>& arrows,int i, int sigma_i) const {
		PartialAutomorphism result=*this;
		result.assign(i,sigma_i);
	 	for (auto i=arrows.begin();i!=arrows.end();) {
	 		switch (result.match_arrow(*i,arrows)) {
	 			case AssignResult::CHANGED : i=arrows.begin(); break;
	 			case AssignResult::INCOMPATIBLE: return nullopt;
	 			case AssignResult::NOCHANGE : ++i;
	 		}
		}
		return result;
	}
	int first_unassigned() const {
		return find(images.begin(),images.end(),UNASSIGNED)-images.begin();
	}
	
	bool complete_permutation() const {
		return find(images.begin(),images.end(),UNASSIGNED)==images.end();
	}


};

template<typename Arrow>
list<vector<int>> automorphisms(const TreeBase<Arrow>& tree) {
	PartialAutomorphism partial{tree.node_hash()};
	return partial.enlargements(tree);
}

#endif	
