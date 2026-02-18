#pragma once

#include <RcppEigen.h>
#include <vinecopulib.hpp>
#include <vinecopulib-wrappers.hpp>
#include <rvine_trees/rvine_trees.hpp>

namespace vinecopulib
{

  inline Rcpp::List rvine_trees_wrap(const RVineTrees &rvine_trees)
  {
    auto [order, struct_array] = rvine_trees.to_struct_array();
    auto struct_array_r = struct_array_wrap(struct_array);

    Rcpp::List rvt = Rcpp::List::create(
        Rcpp::Named("order") = order,
        Rcpp::Named("struct_array") = struct_array_r,
        Rcpp::Named("d") = rvine_trees.get_dim(),
        Rcpp::Named("trunc_lvl") = rvine_trees.get_trunc_lvl());

    rvt.attr("class") = Rcpp::CharacterVector{"rvine_structure", "list"};
    return rvt;
  }

}
