#pragma once
#include <algorithm>
#include <iomanip>
#include <set>
#include <map>
#include <vector>
#include <vinecopulib/misc/tools_stl.hpp>

namespace vinecopulib
{

  namespace tools_stl
  {

    template <class T>
    std::set<T>
    intersect(const std::set<T> &x, const std::set<T> &y)
    {
      std::set<T> out;
      std::set_intersection(
          x.begin(), x.end(), y.begin(), y.end(), std::inserter(out, out.begin()));
      return out;
    }

    template <class T>
    std::set<T>
    set_diff(const std::set<T> &x, const std::set<T> &y)
    {
      std::set<T> out;
      std::set_difference(
          x.begin(), x.end(), y.begin(), y.end(), std::inserter(out, out.begin()));
      return out;
    }

    template <typename T>
    class UnionFind
    {
    public:
      void add(const T &u) { parent_[u] = u; }

      T find(T u)
      {
        while (parent_[u] != u)
        {
          parent_[u] = parent_[parent_[u]]; // path compression
          u = parent_[u];
        }
        return u;
      }

      bool unite(T u, T v)
      {
        T pu = find(u);
        T pv = find(v);
        if (pu == pv)
          return false;
        parent_[pu] = pv;
        return true;
      }

    private:
      std::map<T, T> parent_;
    };
  }
}
