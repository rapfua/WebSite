#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Part 1: Add Edges Based on Gaussian Weight Matrix
// [[Rcpp::export]]
Rcpp::DataFrame add_edges(const arma::mat& W) {
  int n = W.n_rows;
  std::vector<int> from;
  std::vector<int> to;
  std::vector<double> weights;
  
  for (int i = 0; i < n; ++i) {
    if (i % 10 == 0) {
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