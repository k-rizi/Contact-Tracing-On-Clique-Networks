#include <Rcpp.h>
using namespace Rcpp;

List independent_cascade_pre_quarantine_cpp(NumericMatrix adj_matrix, double p, int seed_node, double alpha, int n_iter = 100) {
  int n = adj_matrix.nrow();
  std::set<int> activated_nodes {seed_node};
  std::set<int> new_activated_nodes {seed_node};
  std::set<int> quarantined_nodes;
  
  for (int iter = 0; iter < n_iter; ++iter) {
    if (new_activated_nodes.empty()) {
      break;
    }
    
    std::set<int> next_new_activated_nodes;
    
    for (int node : new_activated_nodes) {
      for (int neighbor = 0; neighbor < n; ++neighbor) {
        if (adj_matrix(node, neighbor) != 0 && activated_nodes.count(neighbor) == 0) {
          // Quarantine the neighbor with probability alpha
          if (R::runif(0, 1) < alpha) {
            quarantined_nodes.insert(neighbor);
          }
          // If the neighbor is not quarantined, attempt to activate it
          else if (quarantined_nodes.count(neighbor) == 0 && R::runif(0, 1) < p) {
            next_new_activated_nodes.insert(neighbor);
          }
        }
      }
    }
    
    activated_nodes.insert(next_new_activated_nodes.begin(), next_new_activated_nodes.end());
    new_activated_nodes.swap(next_new_activated_nodes);
    next_new_activated_nodes.clear();
  }
  
  return List::create(Named("num_active") = activated_nodes.size(),
                      Named("num_quarantine") = quarantined_nodes.size());
}

// [[Rcpp::export]]
DataFrame run_simulations(NumericMatrix adj_matrix, double p, double alpha, int M, int n_iter = 100) {
  int n = adj_matrix.nrow();
  IntegerVector num_active(M);
  IntegerVector num_quarantine(M);
  
  for (int i = 0; i < M; ++i) {
    // Choose a random seed node
    int seed_node = R::runif(0, 1) * n;
    
    // Run the simulation
    List result = independent_cascade_pre_quarantine_cpp(adj_matrix, p, seed_node, alpha, n_iter);
    
    // Store the number of active and quarantine nodes
    num_active[i] = result["num_active"];
    num_quarantine[i] = result["num_quarantine"];
  }
  
  return DataFrame::create(Named("num_active") = num_active,
                           Named("num_quarantine") = num_quarantine);
}