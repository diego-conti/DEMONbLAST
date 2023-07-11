#include "diagramprocessor.h"

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
DiagramProcessor::DiagramProcessor(lie_algebra_list_tag) : processor{new DiagramProcessorListOfLieAlgebras()} {} 

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

