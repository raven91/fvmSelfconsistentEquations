//
// Created by Nikita Kruk on 2019-06-03.
//

#ifndef FVMSELFCONSISTENTEQUATIONS_SOLVERHOMOGENEOUSZEROLAG_HPP
#define FVMSELFCONSISTENTEQUATIONS_SOLVERHOMOGENEOUSZEROLAG_HPP

#include "Definitions.hpp"
#include "HashTuple.hpp"
#include <vector>
#include <unordered_map>
#include <tuple>

class SolverHomogeneousZeroLag
{
 public:

  SolverHomogeneousZeroLag(int n_sigma = 100, int n_D_phi = 100);
  ~SolverHomogeneousZeroLag();

  void FindOrderParameter();

 private:

  int n_sigma_;
  int n_D_phi_;
  Real sigma_min_;
  Real sigma_max_;
  Real dsigma_;
  Real D_phi_min_;
  Real D_phi_max_;
  Real dD_phi_;
  std::unordered_map<std::tuple<Real, Real>, Real, hash_tuple::Hash<std::tuple<Real, Real>>>
      order_parameters_; // (\sigma,D_\varphi) -> R

};

#endif //FVMSELFCONSISTENTEQUATIONS_SOLVERHOMOGENEOUSZEROLAG_HPP
