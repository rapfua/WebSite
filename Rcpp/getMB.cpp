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
  std::unordered_map<std::pair<int, int>, double, pair_hash> edgesMB;
  
  // Constructor: Create a graph from a DataFrame of edges
  graph(Rcpp::DataFrame df,bool approx = false) {
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
  
  // Displays the graph structure
  void display_graph() {
    for (int i = 0; i < n; i++) {
      std::cout << "Node " << i << " has edges: ";
      for (const auto& e : nodes[i].edges) {
        std::cout << "(" << e.u << " -> " << e.v << ", weight: " << e.w << ") ";
      }
      std::cout << std::endl;
    }
  }
  
  // Builds the metric backbone
  void getMB(bool approx = false) {
    std::vector<int> nodesToTry(n);
    std::iota(nodesToTry.begin(), nodesToTry.end(), 0);
    
    if (approx) {
      int cnt = 2 * std::log(n) + 1;
      std::random_device rd;
      std::default_random_engine rng(rd());
      std::shuffle(nodesToTry.begin(), nodesToTry.end(), rng);
      nodesToTry.resize(cnt);
    }
    
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
void create_and_display_graph_from_df(Rcpp::DataFrame df) {
  graph G(df);
  G.display_graph();
}




// [[Rcpp::export]]
arma::mat get_metric_backbone_adjacency_matrix(Rcpp::DataFrame df, bool approx = false) {
  // Create a graph instance from the DataFrame
  graph G(df, approx);
  
  // Return the adjacency matrix of the metric backbone
  return G.get_MB_adjacency_matrix();
}




// [[Rcpp::export]]
arma::mat create_adjacency_matrix(int n, Rcpp::DataFrame edgesMB) {
  // n is the number of nodes in the graph (MB)
  
  // Extracting the edges from the DataFrame
  Rcpp::IntegerVector from = edgesMB["from"];
  Rcpp::IntegerVector to = edgesMB["to"];
  Rcpp::NumericVector weight = edgesMB["weight"];
  
  // Initialize the adjacency matrix with zeros
  arma::mat A = arma::zeros<arma::mat>(n, n);
  
  // Fill in the adjacency matrix based on the edge list
  for (int i = 0; i < from.size(); ++i) {
    int u = from[i];
    int v = to[i];
    double w = weight[i];
    
    // Set the adjacency matrix entries (assuming undirected graph)
    A(u, v) = w;
    A(v, u) = w;  // Since the graph is undirected
  }
  
  return A;
}