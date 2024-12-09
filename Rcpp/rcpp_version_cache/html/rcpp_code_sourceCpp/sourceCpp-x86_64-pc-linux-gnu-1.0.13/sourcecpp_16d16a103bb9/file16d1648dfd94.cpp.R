`.sourceCpp_1_DLLInfo` <- dyn.load('/cloud/project/Rcpp/rcpp_version_cache/html/rcpp_code_sourceCpp/sourceCpp-x86_64-pc-linux-gnu-1.0.13/sourcecpp_16d16a103bb9/sourceCpp_2.so')

get_metric_backbone_adjacency_matrix <- Rcpp:::sourceCppFunction(function(df) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_get_metric_backbone_adjacency_matrix')
add_edges <- Rcpp:::sourceCppFunction(function(W) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_add_edges')
compute_ARI <- Rcpp:::sourceCppFunction(function(true_labels, pred_labels) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_compute_ARI')

rm(`.sourceCpp_1_DLLInfo`)
