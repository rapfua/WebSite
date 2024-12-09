#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <queue>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <random>
#include <fstream>
#include <cmath>
#include <cassert>

using namespace Rcpp;

// Structure to represent an edge in the graph
struct edge {
  int u, v; 
  double w; 
  edge() {}
  edge(int _u, int _v, double _w) : u(_u), v(_v), w(_w) {}
};

// Structure to represent a node containing edges
struct node {
  std::vector<edge> edges;
};

// Helper for hashing pairs
struct pair_hash {
  inline std::size_t operator()(const std::pair<int,int>& v) const {
    return v.first * 31 + v.second;
  }
};

// Graph structure
class graph {
public:
  int n, m;
  std::vector<node> nodes;
  
  //std::unordered_map<int, int, double, pair_hash> edgesMB;
  std::unordered_map<std::pair<int, int>, double, pair_hash> edgesMB;

  
  // Constructor: Create a graph from a DataFrame of edges
  graph(Rcpp::DataFrame df) {
    Rcpp::IntegerVector from = df["from"];
    Rcpp::IntegerVector to = df["to"];
    Rcpp::NumericVector weight = df["weight"];
    

     
    m = from.size();
    n = std::max(*std::max_element(from.begin(), from.end()), *std::max_element(to.begin(), to.end())) + 1;  // Number of nodes
    nodes.resize(n);
    
    for (int i = 0; i < m; i++) {
      int u = from[i];
      int v = to[i];
      double w = weight[i];
      add_edge(u, v, w);
    }


    getMB();
  }
  
  // Returns the adjacency matrix of the metric backbone
  arma::mat get_MB_adjacency_matrix() {
    // Initialize an n x n matrix with zeros
    arma::mat A = arma::zeros<arma::mat>(n, n);
    
    // Fill the adjacency matrix using the edgesMB map
    for (const auto& [edge_pair, weight] : edgesMB) {

      // std::cout << weight << std::endl;
      int u = edge_pair.first;
      int v = edge_pair.second;
      A(u, v) = weight;
      A(v, u) = weight;  // Since the graph is undirected
    }
    
    return A;
  }
  
  
  // Adds an edge between two nodes
  void add_edge(int u, int v, double w) {
    nodes[u].edges.emplace_back(u, v, w);
		nodes[v].edges.emplace_back(v, u, w);

  }
  
  
  // Builds the metric backbone
  void getMB() {
    std::vector<int> nodesToTry(n);
    std::iota(nodesToTry.begin(), nodesToTry.end(), 0);
    
    for (int s : nodesToTry) {
      std::vector<std::pair<int, double>> parent(n, {-1, -1});
      getSPT(s, parent);
      for (int u = 0; u < n; ++u) {
        auto [v, w] = parent[u];
        if (v != -1) {
          edgesMB[{u, v}] = w;
        }
      }
    }
  }
  
  // Builds a shortest path tree (SPT) using Dijkstra's algorithm
  void getSPT(int s, std::vector<std::pair<int, double>>& parent) {
    std::vector<double> dist(n, 1e15);
    std::vector<bool> visited(n, false);
    std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<>> pq;
    dist[s] = 0.0;
    pq.emplace(0.0, s);
    
    while (!pq.empty()) {
      int cur = pq.top().second;
      pq.pop();
      if (visited[cur]) continue;
      visited[cur] = true;
      
      for (const auto& e : nodes[cur].edges) {
        if (dist[e.v] > dist[cur] + e.w) {
          dist[e.v] = dist[cur] + e.w;
          parent[e.v] = {e.u, e.w};
          pq.emplace(dist[e.v], e.v);
        }
      }
    }
  }
  

};




// [[Rcpp::export]]
arma::mat get_metric_backbone_adjacency_matrix(Rcpp::DataFrame df) {
  // Create a graph instance from the DataFrame
  graph G(df);
  
  // Return the adjacency matrix of the metric backbone
  return G.get_MB_adjacency_matrix();
}



// Part 1: Add Edges Based on Gaussian Weight Matrix
// [[Rcpp::export]]
Rcpp::DataFrame add_edges(const arma::mat& W) {
  int n = W.n_rows;
  std::vector<int> from;
  std::vector<int> to;
  std::vector<double> weights;
  
  for (int i = 0; i < n; ++i) {
    if (i % 100 == 0) {
      Rcpp::Rcout << "Processing node " << i << std::endl;
    }
    
    for (int j = i + 1; j < n; ++j) {
      double w = W(i, j);
      if (w > 0) {
        from.push_back(i);
        to.push_back(j);
        weights.push_back(1.0 / w - 1.0);
      }
      
    }
    
  }
  
  return Rcpp::DataFrame::create(Rcpp::Named("from") = from,
                                 Rcpp::Named("to") = to,
                                 Rcpp::Named("weight") = weights);
}


// Part 2: Compute Adjusted Rand Index (ARI)
// [[Rcpp::export]]
double compute_ARI(Rcpp::IntegerVector true_labels, Rcpp::IntegerVector pred_labels) {
  Rcpp::Environment mclust("package:mclust");
  Rcpp::Function adjustedRandIndex = mclust["adjustedRandIndex"];
  
  // Call R's adjustedRandIndex function
  Rcpp::NumericVector result = adjustedRandIndex(true_labels, pred_labels);
  
  // Return the result as a double
  return Rcpp::as<double>(result);
}



#include <Rcpp.h>
#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// get_metric_backbone_adjacency_matrix
arma::mat get_metric_backbone_adjacency_matrix(Rcpp::DataFrame df);
RcppExport SEXP sourceCpp_1_get_metric_backbone_adjacency_matrix(SEXP dfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type df(dfSEXP);
    rcpp_result_gen = Rcpp::wrap(get_metric_backbone_adjacency_matrix(df));
    return rcpp_result_gen;
END_RCPP
}
// add_edges
Rcpp::DataFrame add_edges(const arma::mat& W);
RcppExport SEXP sourceCpp_1_add_edges(SEXP WSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type W(WSEXP);
    rcpp_result_gen = Rcpp::wrap(add_edges(W));
    return rcpp_result_gen;
END_RCPP
}
// compute_ARI
double compute_ARI(Rcpp::IntegerVector true_labels, Rcpp::IntegerVector pred_labels);
RcppExport SEXP sourceCpp_1_compute_ARI(SEXP true_labelsSEXP, SEXP pred_labelsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type true_labels(true_labelsSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type pred_labels(pred_labelsSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_ARI(true_labels, pred_labels));
    return rcpp_result_gen;
END_RCPP
}
