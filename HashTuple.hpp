//
// Created by Nikita Kruk on 2019-06-03.
//

#ifndef FVMSELFCONSISTENTEQUATIONS_HASHTUPLE_HPP
#define FVMSELFCONSISTENTEQUATIONS_HASHTUPLE_HPP

#include <cstddef>
#include <functional>
#include <tuple>

namespace hash_tuple
{
  template<typename TT>
  struct Hash
  {
    std::size_t operator()(TT const &tt) const
    {
      return std::hash<TT>()(tt);
    }
  };

  namespace
  {
    template<class T>
    inline void HashCombine(std::size_t &seed, T const &v)
    {
      seed ^= hash_tuple::Hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
  }

  namespace
  {
    template<class Tuple, std::size_t Index = std::tuple_size<Tuple>::value - 1>
    struct HashValueImpl
    {
      static void Apply(std::size_t &seed, Tuple const &tuple)
      {
        HashValueImpl<Tuple, Index - 1>::Apply(seed, tuple);
        HashCombine(seed, std::get<Index>(tuple));
      }
    };

    template<class Tuple>
    struct HashValueImpl<Tuple, 0>
    {
      static void Apply(std::size_t &seed, Tuple const &tuple)
      {
        HashCombine(seed, std::get<0>(tuple));
      }
    };
  }

  template<typename ... TT>
  struct Hash<std::tuple<TT...>>
  {
    std::size_t operator()(std::tuple<TT...> const &tt) const
    {
      std::size_t seed = 0;
      HashValueImpl<std::tuple<TT...>>::Apply(seed, tt);
      return seed;
    }
  };
}

#endif //FVMSELFCONSISTENTEQUATIONS_HASHTUPLE_HPP
