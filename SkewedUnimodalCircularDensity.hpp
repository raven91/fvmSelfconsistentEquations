//
// Created by Nikita Kruk on 01.01.20.
//

#ifndef FVMSELFCONSISTENTEQUATIONS_SKEWEDUNIMODALCIRCULARDENSITY_HPP
#define FVMSELFCONSISTENTEQUATIONS_SKEWEDUNIMODALCIRCULARDENSITY_HPP

#include "Definitions.hpp"

class SkewedUnimodalCircularDensity
{
 public:

  SkewedUnimodalCircularDensity();
  explicit SkewedUnimodalCircularDensity(Real coupling,
                                         Real lag,
                                         Real diffusion,
                                         Real group_velocity,
                                         Real order_parameter);
  ~SkewedUnimodalCircularDensity();

  void ReinitializeParameters(Real coupling, Real lag, Real diffusion, Real group_velocity, Real order_parameter);
  void IntegrateForRealPart(MultiprecisionReal &real_part) const;
  void IntegrateForImaginaryPart(MultiprecisionReal &imaginary_part) const;

 private:

  Real coupling_;
  Real lag_;
  Real diffusion_;
  Real inverse_diffusion_;
  Real group_velocity_;
  Real order_parameter_;

  int number_of_quadrature_points_;
  std::vector<Real> abstract_quadrature_points_;
  std::vector<Real> abstract_quadrature_weights_;
  int number_of_grid_points_;
  std::vector<Real> decomposed_grid_points_;
  std::vector<MultiprecisionReal> exponent_at_decomposed_point_;
  std::vector<MultiprecisionReal> inverse_exponent_at_decomposed_point_;
  std::vector<MultiprecisionReal> integral_up_to_grid_point_;
  std::vector<MultiprecisionReal> integral_up_to_decomposed_point_;
  MultiprecisionReal normalization_constant_;
  std::vector<Real> cos_of_decomposed_point_;
  std::vector<Real> sin_of_decomposed_point_;

  void E(Real x, MultiprecisionReal &res) const;
  void E(const MultiprecisionReal &x, MultiprecisionReal &res) const;
  void EInv(Real x, MultiprecisionReal &res) const;
  void EInv(const MultiprecisionReal &x, MultiprecisionReal &res) const;
  void IntegrateInverseExponentForGrid(Real lower_limit, Real upper_limit);
  void ComputeDensityFunctionNormalizationConstant();
  void CalculateAuxiliaryExponentsAndIntegrals();

};

#endif //FVMSELFCONSISTENTEQUATIONS_SKEWEDUNIMODALCIRCULARDENSITY_HPP
