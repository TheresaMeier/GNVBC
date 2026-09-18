#pragma once

#include <vinecopulib.hpp>

#include <algorithm>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace gnvbc {

using Edge = vinecopulib::RVineTrees::Edge;
using Tree = vinecopulib::RVineTrees::Tree;

// Mutable state used while completing a forest. It deliberately does not
// retain pointers into RVineTrees' augmented conversion representation.
struct CompletionEdge {
  size_t node1;
  size_t node2;
  Edge edge;
};

class UnionFind {
public:
  explicit UnionFind(size_t n)
    : parent_(n)
    , rank_(n, 0)
  {
    std::iota(parent_.begin(), parent_.end(), 0);
  }

  size_t find(size_t node)
  {
    if (parent_.at(node) != node)
      parent_[node] = find(parent_[node]);
    return parent_[node];
  }

  bool unite(size_t left, size_t right)
  {
    left = find(left);
    right = find(right);
    if (left == right)
      return false;
    if (rank_[left] < rank_[right])
      std::swap(left, right);
    parent_[right] = left;
    if (rank_[left] == rank_[right])
      ++rank_[left];
    return true;
  }

private:
  std::vector<size_t> parent_;
  std::vector<size_t> rank_;
};

inline std::vector<size_t>
insert_sorted(std::vector<size_t> values, size_t value)
{
  values.insert(std::lower_bound(values.begin(), values.end(), value), value);
  return values;
}

inline std::string edge_description(const Edge& edge)
{
  std::string result = "(" + std::to_string(edge.a) + ", " +
                       std::to_string(edge.b) + ";";
  for (auto value : edge.C)
    result += " " + std::to_string(value);
  return result + ")";
}

inline void validate_edge(const Edge& edge, size_t tree, size_t d)
{
  if (edge.a == 0 || edge.a > d || edge.b == 0 || edge.b > d ||
      edge.a == edge.b || edge.C.size() != tree ||
      !std::is_sorted(edge.C.begin(), edge.C.end()) ||
      std::adjacent_find(edge.C.begin(), edge.C.end()) != edge.C.end() ||
      std::find(edge.C.begin(), edge.C.end(), edge.a) != edge.C.end() ||
      std::find(edge.C.begin(), edge.C.end(), edge.b) != edge.C.end()) {
    throw std::runtime_error("Invalid inherited edge in tree " +
                             std::to_string(tree) + ": " +
                             edge_description(edge));
  }
}

inline std::map<std::pair<size_t, std::vector<size_t>>, size_t>
make_node_lookup(const std::vector<CompletionEdge>& previous)
{
  std::map<std::pair<size_t, std::vector<size_t>>, size_t> lookup;
  for (size_t node = 0; node < previous.size(); ++node) {
    const Edge& edge = previous[node].edge;
    lookup[{ edge.a, insert_sorted(edge.C, edge.b) }] = node;
    lookup[{ edge.b, insert_sorted(edge.C, edge.a) }] = node;
  }
  return lookup;
}

inline std::vector<size_t> intersection(const std::vector<size_t>& left,
                                        const std::vector<size_t>& right)
{
  std::vector<size_t> result;
  std::set_intersection(left.begin(), left.end(), right.begin(), right.end(),
                        std::back_inserter(result));
  return result;
}

inline std::vector<size_t> difference(const std::vector<size_t>& left,
                                      const std::vector<size_t>& right)
{
  std::vector<size_t> result;
  std::set_difference(left.begin(), left.end(), right.begin(), right.end(),
                      std::back_inserter(result));
  return result;
}

inline std::vector<CompletionEdge>
complete_tree(size_t tree, size_t d, const Tree& inherited,
              const std::vector<CompletionEdge>& previous)
{
  const size_t target_edges = d - 1 - tree;
  if (inherited.size() > target_edges) {
    throw std::runtime_error("Tree " + std::to_string(tree) + " has " +
                             std::to_string(inherited.size()) +
                             " inherited edges; expected at most " +
                             std::to_string(target_edges) + ".");
  }

  std::vector<std::vector<size_t>> node_variables;
  std::vector<CompletionEdge> result;
  if (tree == 0) {
    node_variables.resize(d);
    for (size_t node = 0; node < d; ++node)
      node_variables[node] = { node + 1 };
    for (const Edge& edge : inherited) {
      validate_edge(edge, tree, d);
      result.push_back({ edge.a - 1, edge.b - 1, edge });
    }
  } else {
    if (previous.size() != d - tree) {
      throw std::runtime_error("Cannot complete tree " + std::to_string(tree) +
                               ": tree below is incomplete.");
    }
    node_variables.reserve(previous.size());
    for (const auto& node : previous)
      node_variables.push_back(node.edge.all_indices);
    const auto lookup = make_node_lookup(previous);
    for (const Edge& edge : inherited) {
      validate_edge(edge, tree, d);
      const auto left = lookup.find({ edge.a, edge.C });
      const auto right = lookup.find({ edge.b, edge.C });
      if (left == lookup.end() || right == lookup.end() ||
          left->second == right->second) {
        throw std::runtime_error("Proximity condition violated by inherited "
                                 "edge in tree " + std::to_string(tree) +
                                 ": " + edge_description(edge));
      }
      result.push_back({ left->second, right->second, edge });
    }
  }

  UnionFind components(node_variables.size());
  std::set<std::pair<size_t, size_t>> connected;
  for (const auto& edge : result) {
    const auto key = std::minmax(edge.node1, edge.node2);
    if (!connected.insert(key).second || !components.unite(edge.node1, edge.node2)) {
      throw std::runtime_error("Inherited forest contains a cycle in tree " +
                               std::to_string(tree) + ".");
    }
  }

  std::vector<CompletionEdge> candidates;
  for (size_t left = 0; left < node_variables.size(); ++left) {
    for (size_t right = left + 1; right < node_variables.size(); ++right) {
      if (connected.count({ left, right }) != 0)
        continue;
      const auto common = intersection(node_variables[left], node_variables[right]);
      if (common.size() != node_variables[left].size() - 1 ||
          common.size() != node_variables[right].size() - 1)
        continue;
      const auto left_only = difference(node_variables[left], common);
      const auto right_only = difference(node_variables[right], common);
      if (left_only.size() == 1 && right_only.size() == 1)
        candidates.push_back(
          { left, right, Edge(left_only[0], right_only[0], common) });
    }
  }

  // Match the former completion algorithm's deterministic reverse candidate
  // traversal while avoiding upstream conversion internals.
  while (result.size() < target_edges && !candidates.empty()) {
    const CompletionEdge candidate = candidates.back();
    candidates.pop_back();
    if (components.unite(candidate.node1, candidate.node2))
      result.push_back(candidate);
  }
  if (result.size() != target_edges) {
    throw std::runtime_error("Could not complete tree " + std::to_string(tree) +
                             " from the inherited forest.");
  }
  return result;
}

inline std::vector<Tree>
complete_trees(size_t d, const std::vector<Tree>& inherited)
{
  if (d < 2 || inherited.empty() || inherited.size() >= d)
    throw std::runtime_error("A non-degenerate, truncated vine is required for merging.");

  std::vector<Tree> completed(inherited.size());
  std::vector<CompletionEdge> previous;
  for (size_t tree = 0; tree < inherited.size(); ++tree) {
    previous = complete_tree(tree, d, inherited[tree], previous);
    completed[tree].reserve(previous.size());
    for (const auto& edge : previous)
      completed[tree].push_back(edge.edge);
  }
  return completed;
}

inline Edge remap_edge(const Edge& edge, const std::vector<size_t>& map)
{
  const auto remap = [&map](size_t variable) {
    if (variable == 0 || variable > map.size() || map[variable - 1] == 0)
      throw std::runtime_error("Local-to-global map does not cover variable " +
                               std::to_string(variable) + ".");
    return map[variable - 1];
  };
  std::vector<size_t> conditioning;
  conditioning.reserve(edge.C.size());
  for (size_t variable : edge.C)
    conditioning.push_back(remap(variable));
  return Edge(remap(edge.a), remap(edge.b), std::move(conditioning));
}

inline vinecopulib::RVineStructure
merge_components(size_t d, const std::vector<vinecopulib::RVineTrees>& components,
                 const std::vector<std::vector<size_t>>& maps)
{
  if (d < 2 || components.empty() || components.size() != maps.size())
    throw std::runtime_error("At least one non-degenerate component and map are required.");

  size_t trunc_lvl = 0;
  std::vector<bool> global_seen(d + 1, false);
  for (size_t component = 0; component < components.size(); ++component) {
    if (components[component].get_dim() < 2 ||
        components[component].get_trunc_lvl() == 0 ||
        maps[component].size() != components[component].get_dim()) {
      throw std::runtime_error("Each component must be non-degenerate and have "
                               "one global label per local variable.");
    }
    std::set<size_t> local_global;
    for (size_t label : maps[component]) {
      if (label == 0 || label > d || !local_global.insert(label).second)
        throw std::runtime_error("Invalid local-to-global variable map.");
      global_seen[label] = true;
    }
    trunc_lvl = std::max(trunc_lvl, components[component].get_trunc_lvl());
  }
  for (size_t variable = 1; variable <= d; ++variable) {
    if (!global_seen[variable])
      throw std::runtime_error("Local-to-global maps omit global variable " +
                               std::to_string(variable) + ".");
  }

  std::vector<Tree> inherited(trunc_lvl);
  for (size_t component = 0; component < components.size(); ++component) {
    const auto& trees = components[component].get_trees();
    for (size_t tree = 0; tree < trees.size(); ++tree) {
      for (const auto& edge : trees[tree])
        inherited[tree].push_back(remap_edge(edge, maps[component]));
    }
  }

  const auto completed = complete_trees(d, inherited);
  return vinecopulib::RVineStructure(vinecopulib::RVineTrees(d, completed));
}

} // namespace gnvbc
