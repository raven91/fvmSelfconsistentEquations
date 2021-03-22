//
// Created by Nikita Kruk on 2019-06-03.
//

#ifndef FVMSELFCONSISTENTEQUATIONS_SOLVERHOMOGENEOUSNONZEROLAG_HPP
#define FVMSELFCONSISTENTEQUATIONS_SOLVERHOMOGENEOUSNONZEROLAG_HPP

#include "Definitions.hpp"
#include "HashTuple.hpp"
#include "ThreadForParameterSpan.hpp"

#include <unordered_map>
#include <tuple>
#include <vector>
#include <random>
#include <fstream>
#include <array>

class SolverHomogeneousNonzeroLag
{
 public:

  SolverHomogeneousNonzeroLag(Thread *thread);
  ~SolverHomogeneousNonzeroLag();

//  void FindOrderParameter_partitive();
  void FindOrderParameter_coupled();
  void FindModeOfDensityFunction();
  void ComputeSkewnessOfDensityFunction();

 private:

  Thread *thread_;

  int number_of_phase_grid_points_;
  int n_sigma_;
  int n_D_phi_;
  int n_alpha_;
  Real sigma_min_;
  Real sigma_max_;
  Real dsigma_;
  Real D_phi_min_;
  Real D_phi_max_;
  Real dD_phi_;
  Real alpha_min_;
  Real alpha_max_;
  Real dalpha_;
  std::unordered_map<std::tuple<Real, Real, Real>, Real, hash_tuple::Hash<std::tuple<Real, Real, Real>>>
      order_parameters_; // (\sigma,D_\varphi,\alpha) -> R
  std::unordered_map<std::tuple<Real, Real, Real>, Real, hash_tuple::Hash<std::tuple<Real, Real, Real>>>
      velocities_; // (\sigma,D_\varphi,\alpha) -> v

  std::vector<MultiprecisionReal> exponent_at_current_point_;
  std::vector<MultiprecisionReal> inverse_exponent_at_current_point_;
  std::vector<MultiprecisionReal> integral_up_to_current_point_;
  std::mt19937 mersenne_twister_generator_;
  std::uniform_real_distribution<Real> uniform_real_distribution_;

  void CreateParameterFiles(std::ifstream &input_parameter_file, std::ofstream &output_parameter_file) const;
  void CreateOutputFiles(std::ofstream &order_parameter_file, std::ofstream &velocity_file) const;
  void CreateInputFiles(std::ifstream &order_parameter_file, std::ifstream &velocity_file) const;
  void CreateOutputFile(std::ofstream &file, const std::string &subfolder, const std::string &base_name) const;
  void DefineParametersForComparisonToFvm(std::vector<std::array<Real, 3>> &comparison_to_fvm_parameters) const;

  // complementary functions inside the solution
//  Real E(Real phi, Real v, Real R, Real D_phi, Real sigma, Real alpha) const;
  void E(Real phi, Real v, Real R, Real D_phi, Real sigma, Real alpha, MultiprecisionReal &res) const;
//  Real EInv(Real phi, Real v, Real R, Real D_phi, Real sigma, Real alpha) const;
  void EInv(Real phi, Real v, Real R, Real D_phi, Real sigma, Real alpha, MultiprecisionReal &res) const;
//  Real IntegrateInverseExponentUpTo(Real phi, Real v, Real R, Real D_phi, Real sigma, Real alpha) const;
//  void IntegrateInverseExponentOverRange(Real lower_limit,
//                                         Real upper_limit,
//                                         Real v,
//                                         Real R,
//                                         Real D_phi,
//                                         Real sigma,
//                                         Real alpha);
  void IntegrateInverseExponentOverRange(Real lower_limit, Real upper_limit);
//  Real ComputeDensityFunctionNormalizationConstant(Real v, Real R, Real D_phi, Real sigma, Real alpha) const;
  void DensityFunctionNormalization(Real v,
                                    Real D_phi,
                                    MultiprecisionReal &normalization) const;
//  Real ImaginaryPart(Real v, Real R, Real D_phi, Real sigma, Real alpha, Real normalization) const;
  void ImaginaryPart(Real v,
                     Real D_phi,
                     const MultiprecisionReal &normalization,
                     MultiprecisionReal &imaginary_part) const;
//  Real RealPart(Real v, Real R, Real D_phi, Real sigma, Real alpha, Real normalization) const;
  void RealPart(Real v,
                Real D_phi,
                const MultiprecisionReal &normalization,
                MultiprecisionReal &real_part) const;
  void CalculateAuxiliaryExponentsAndIntegrals(Real v, Real R, Real D_phi, Real sigma, Real alpha);
  bool IsInRange(Real p, Real p_min, Real p_max) const;

};

#endif //FVMSELFCONSISTENTEQUATIONS_SOLVERHOMOGENEOUSNONZEROLAG_HPP
