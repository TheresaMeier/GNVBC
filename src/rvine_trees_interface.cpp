#include "rvine_trees_wrappers.hpp"

using namespace vinecopulib;

// [[Rcpp::export]]
Rcpp::List merge_rvine_structures(Rcpp::List rvine_structure_list)
{
  // Check input length
  if (rvine_structure_list.size() == 0)
  {
    Rcpp::stop("At least one rvine_structure must be provided.");
  }

  // Convert to vector of RVineTrees
  std::vector<RVineTrees> vines;
  for (int i = 0; i < rvine_structure_list.size(); ++i)
  {
    const Rcpp::List &r_struct = rvine_structure_list[i];
    size_t trunc_lvl = r_struct["trunc_lvl"];
    std::vector<size_t> order = r_struct["order"];
    TriangularArray<size_t> struct_array = struct_array_wrap(r_struct["struct_array"], trunc_lvl);
    vines.push_back(RVineTrees(order, struct_array));
  }

  // Merge vines
  RVineTrees merged_vine(vines);

  // Convert to struct_array (with fill_missing = true)
  auto [merged_order, merged_array] = merged_vine.to_struct_array(true);

  // Create a valid RVineStructure
  RVineStructure merged_structure(merged_order, merged_array);

  // Return wrapped R object
  return rvine_structure_wrap(merged_structure);
}
