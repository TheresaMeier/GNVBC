#include <RcppEigen.h>
#include <vinecopulib.hpp>
#include <vinecopulib-wrappers.hpp>

#include <gnvbc/rvine_merge.hpp>

using namespace vinecopulib;

// [[Rcpp::export]]
Rcpp::List merge_rvine_structures(Rcpp::List rvine_structure_list,
                                  Rcpp::List local_to_global_maps,
                                  int global_dim)
{
  if (rvine_structure_list.size() == 0 ||
      rvine_structure_list.size() != local_to_global_maps.size()) {
    Rcpp::stop("At least one structure and one map per structure are required.");
  }
  if (global_dim < 2)
    Rcpp::stop("The merged vine dimension must be at least two.");

  std::vector<RVineTrees> components;
  std::vector<std::vector<size_t>> maps;
  components.reserve(rvine_structure_list.size());
  maps.reserve(rvine_structure_list.size());
  for (int i = 0; i < rvine_structure_list.size(); ++i) {
    const Rcpp::List& r_struct = rvine_structure_list[i];
    RVineStructure structure = rvine_structure_wrap(r_struct, true);
    components.push_back(structure.get_trees());
    maps.push_back(Rcpp::as<std::vector<size_t>>(local_to_global_maps[i]));
  }

  RVineStructure merged_structure =
    gnvbc::merge_components(static_cast<size_t>(global_dim), components, maps);
  return rvine_structure_wrap(merged_structure);
}
