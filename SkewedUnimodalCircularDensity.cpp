//
// Created by Nikita Kruk on 01.01.20.
//

#include "SkewedUnimodalCircularDensity.hpp"

#define TRAPEZOIDAL_RULE
//#define GAUSS_LEGENDRE_QUADRATURE

#include <cmath>
#include <iostream>
#include <numeric> // std::inner_product

SkewedUnimodalCircularDensity::SkewedUnimodalCircularDensity() :
    coupling_(0.0),
    lag_(0.0),
    diffusion_(0.0),
    inverse_diffusion_(0.0),
    group_velocity_(0.0),
    order_parameter_(0.0),
    number_of_quadrature_points_(1),
    abstract_quadrature_points_(number_of_quadrature_points_, 0.0),
    abstract_quadrature_weights_(number_of_quadrature_points_, 0.0),
    number_of_grid_points_(10000),
    decomposed_grid_points_(number_of_grid_points_ * number_of_quadrature_points_, 0.0),
    exponent_at_decomposed_point_(number_of_grid_points_ * number_of_quadrature_points_, MultiprecisionReal(0.0)),
    inverse_exponent_at_decomposed_point_(number_of_grid_points_ * number_of_quadrature_points_,
                                          MultiprecisionReal(0.0)),
    integral_up_to_grid_point_(number_of_grid_points_ + 1, MultiprecisionReal(0.0)),
    integral_up_to_decomposed_point_(number_of_grid_points_ * number_of_quadrature_points_, MultiprecisionReal(0.0)),
    normalization_constant_(0.0),
    cos_of_decomposed_point_(number_of_grid_points_ * number_of_quadrature_points_, 0.0),
    sin_of_decomposed_point_(number_of_grid_points_ * number_of_quadrature_points_, 0.0)
{
#if defined(TRAPEZOIDAL_RULE)
  const Real dx = kTwoPi / number_of_grid_points_;
  static Real x_cur = 0.0;
  for (int i = 0; i < number_of_grid_points_ + 1; ++i)
  {
    x_cur = i * dx;
    cos_of_decomposed_point_[i] = std::cos(x_cur);
    sin_of_decomposed_point_[i] = std::sin(x_cur);
  } // i
#elif defined(GAUSS_LEGENDRE_QUADRATURE)
  const Real dx = kTwoPi / number_of_grid_points_;
  GenerateGaussLegendrePoints(0.0,
                              dx,
                              number_of_quadrature_points_,
                              abstract_quadrature_points_,
                              abstract_quadrature_weights_);
  for (int i = 0; i < number_of_grid_points_; ++i)
  {
    for (int m = 0; m < number_of_quadrature_points_; ++m)
    {
      decomposed_grid_points_[m + i * number_of_quadrature_points_] = i * dx + abstract_quadrature_points_[m];
      cos_of_decomposed_point_[m + i * number_of_quadrature_points_] =
          std::cos(decomposed_grid_points_[m + i * number_of_quadrature_points_]);
      sin_of_decomposed_point_[m + i * number_of_quadrature_points_] =
          std::sin(decomposed_grid_points_[m + i * number_of_quadrature_points_]);
    } // m
  } // i
#endif
}

SkewedUnimodalCircularDensity::SkewedUnimodalCircularDensity(Real coupling,
                                                             Real lag,
                                                             Real diffusion,
                                                             Real group_velocity,
                                                             Real order_parameter) :
    coupling_(coupling),
    lag_(lag),
    diffusion_(diffusion),
    inverse_diffusion_(1.0 / diffusion),
    group_velocity_(group_velocity),
    order_parameter_(order_parameter),
    number_of_quadrature_points_(3),
    abstract_quadrature_points_(number_of_quadrature_points_, 0.0),
    abstract_quadrature_weights_(number_of_quadrature_points_, 0.0),
    number_of_grid_points_(1000),
    decomposed_grid_points_(number_of_grid_points_ * number_of_quadrature_points_, 0.0),
    exponent_at_decomposed_point_(number_of_grid_points_ * number_of_quadrature_points_, MultiprecisionReal(0.0)),
    inverse_exponent_at_decomposed_point_(number_of_grid_points_ * number_of_quadrature_points_,
                                          MultiprecisionReal(0.0)),
    integral_up_to_grid_point_(number_of_grid_points_ + 1, MultiprecisionReal(0.0)),
    integral_up_to_decomposed_point_(number_of_grid_points_ * number_of_quadrature_points_, MultiprecisionReal(0.0)),
    normalization_constant_(0.0),
    cos_of_decomposed_point_(number_of_grid_points_ * number_of_quadrature_points_, 0.0),
    sin_of_decomposed_point_(number_of_grid_points_ * number_of_quadrature_points_, 0.0)
{
#if defined(TRAPEZOIDAL_RULE)
  const Real dx = kTwoPi / number_of_grid_points_;
  static Real x_cur = 0.0;
  for (int i = 0; i < number_of_grid_points_ + 1; ++i)
  {
    x_cur = i * dx;
    cos_of_decomposed_point_[i] = std::cos(x_cur);
    sin_of_decomposed_point_[i] = std::sin(x_cur);
  } // i
#elif defined(GAUSS_LEGENDRE_QUADRATURE)
  const Real dx = kTwoPi / number_of_grid_points_;
  GenerateGaussLegendrePoints(0.0,
                              dx,
                              number_of_quadrature_points_,
                              abstract_quadrature_points_,
                              abstract_quadrature_weights_);
  for (int i = 0; i < number_of_grid_points_; ++i)
  {
    for (int m = 0; m < number_of_quadrature_points_; ++m)
    {
      decomposed_grid_points_[m + i * number_of_quadrature_points_] = i * dx + abstract_quadrature_points_[m];
      cos_of_decomposed_point_[m + i * number_of_quadrature_points_] =
          std::cos(decomposed_grid_points_[m + i * number_of_quadrature_points_]);
      sin_of_decomposed_point_[m + i * number_of_quadrature_points_] =
          std::sin(decomposed_grid_points_[m + i * number_of_quadrature_points_]);
    } // m
  } // i
#endif

  CalculateAuxiliaryExponentsAndIntegrals();
  ComputeDensityFunctionNormalizationConstant();
}

SkewedUnimodalCircularDensity::~SkewedUnimodalCircularDensity() = default;

void SkewedUnimodalCircularDensity::ReinitializeParameters(Real coupling,
                                                           Real lag,
                                                           Real diffusion,
                                                           Real group_velocity,
                                                           Real order_parameter)
{
  coupling_ = coupling;
  lag_ = lag;
  diffusion_ = diffusion;
  inverse_diffusion_ = 1.0 / diffusion;
  group_velocity_ = group_velocity;
  order_parameter_ = order_parameter;

  /*std::fill(exponent_at_decomposed_point_.begin(), exponent_at_decomposed_point_.end(), MultiprecisionReal(0.0));
  std::fill(inverse_exponent_at_decomposed_point_.begin(),
            inverse_exponent_at_decomposed_point_.end(),
            MultiprecisionReal(0.0));
  std::fill(integral_up_to_grid_point_.begin(), integral_up_to_grid_point_.end(), MultiprecisionReal(0.0));
  std::fill(integral_up_to_decomposed_point_.begin(), integral_up_to_decomposed_point_.end(), MultiprecisionReal(0.0));
  normalization_constant_ = 0.0;*/

  CalculateAuxiliaryExponentsAndIntegrals();
  ComputeDensityFunctionNormalizationConstant();
}

void SkewedUnimodalCircularDensity::CalculateAuxiliaryExponentsAndIntegrals()
{
#if defined(TRAPEZOIDAL_RULE)
  const Real dx = kTwoPi / number_of_grid_points_;
  static Real x_cur = 0.0;
  for (int i = 0; i < number_of_grid_points_ + 1; ++i)
  {
    x_cur = i * dx;
    E(x_cur, exponent_at_decomposed_point_[i]);
    EInv(x_cur, inverse_exponent_at_decomposed_point_[i]);
  } // i
#elif defined(GAUSS_LEGENDRE_QUADRATURE)
  for (int i = 0; i < decomposed_grid_points_.size(); ++i)
  {
    E(decomposed_grid_points_[i], exponent_at_decomposed_point_[i]);
    EInv(decomposed_grid_points_[i], inverse_exponent_at_decomposed_point_[i]);
  } // i
#endif
  const Real lower_limit = 0.0, &upper_limit = kTwoPi;
  IntegrateInverseExponentForGrid(lower_limit, upper_limit);
}

void SkewedUnimodalCircularDensity::E(Real x, MultiprecisionReal &res) const
{
  static MultiprecisionReal exp_arg = 0.0;
  exp_arg =
      (-group_velocity_ * x + coupling_ * order_parameter_ * std::cos(x + lag_)) * inverse_diffusion_;
  res = boost::multiprecision::exp(exp_arg);
}

void SkewedUnimodalCircularDensity::E(const MultiprecisionReal &x, MultiprecisionReal &res) const
{
  static MultiprecisionReal exp_arg = 0.0;
  exp_arg = (-group_velocity_ * x
      + coupling_ * order_parameter_ * boost::multiprecision::cos(x + lag_)) * inverse_diffusion_;
  res = boost::multiprecision::exp(exp_arg);
}

void SkewedUnimodalCircularDensity::EInv(Real x, MultiprecisionReal &res) const
{
  static MultiprecisionReal exp_arg = 0.0;
  exp_arg =
      (group_velocity_ * x - coupling_ * order_parameter_ * std::cos(x + lag_)) * inverse_diffusion_;
  res = boost::multiprecision::exp(exp_arg);
}

void SkewedUnimodalCircularDensity::EInv(const MultiprecisionReal &x, MultiprecisionReal &res) const
{
  static MultiprecisionReal exp_arg = 0.0;
  exp_arg = (group_velocity_ * x
      - coupling_ * order_parameter_ * boost::multiprecision::cos(x + lag_)) * inverse_diffusion_;
  res = boost::multiprecision::exp(exp_arg);
}

void SkewedUnimodalCircularDensity::IntegrateInverseExponentForGrid(Real lower_limit, Real upper_limit)
{
  const Real dx = (upper_limit - lower_limit) / number_of_grid_points_;
#if defined(TRAPEZOIDAL_RULE)
  integral_up_to_grid_point_[0] = 0.0;
  for (int i = 0; i < number_of_grid_points_; ++i)
  {
    integral_up_to_grid_point_[i + 1] = integral_up_to_grid_point_[i]
        + (inverse_exponent_at_decomposed_point_[i] + inverse_exponent_at_decomposed_point_[i + 1]) * (0.5 * dx);
  } // i
#elif defined(GAUSS_LEGENDRE_QUADRATURE)
  integral_up_to_grid_point_[0] = 0.0;
  for (int i = 1; i <= number_of_grid_points_; ++i)
  {
    integral_up_to_grid_point_[i] = integral_up_to_grid_point_[i - 1]
        + std::inner_product(abstract_quadrature_weights_.begin(),
                             abstract_quadrature_weights_.end(),
                             &inverse_exponent_at_decomposed_point_[(i - 1) * number_of_quadrature_points_],
                             MultiprecisionReal(0.0));
  } // i
  std::vector<Real> subquadrature_points(number_of_quadrature_points_, 0.0),
      subquadrature_weights(number_of_quadrature_points_, 0.0);
  static MultiprecisionReal inv_exp(0.0);
  for (int i = 0; i < number_of_grid_points_; ++i)
  {
    for (int m = 0; m < number_of_quadrature_points_; ++m)
    {
      GenerateGaussLegendrePoints(i * dx,
                                  decomposed_grid_points_[m + i * number_of_quadrature_points_],
                                  number_of_quadrature_points_,
                                  subquadrature_points,
                                  subquadrature_weights);
      integral_up_to_decomposed_point_[m + i * number_of_quadrature_points_] = integral_up_to_grid_point_[i];
      for (int s = 0; s < subquadrature_points.size(); ++s)
      {
        EInv(subquadrature_points[s], inv_exp);
        integral_up_to_decomposed_point_[m + i * number_of_quadrature_points_] += subquadrature_weights[s] * inv_exp;
      } // s
    } // m
  } // i
#endif
}

void SkewedUnimodalCircularDensity::ComputeDensityFunctionNormalizationConstant()
{
  // auxiliary variables
  const MultiprecisionReal &integral_of_inv_exp_02pi
      = integral_up_to_grid_point_[integral_up_to_grid_point_.size() - 1];
  static MultiprecisionReal periodicity_constraint = 0.0, exp_arg = 0.0;
  exp_arg = kTwoPi * group_velocity_ * inverse_diffusion_;
  periodicity_constraint = boost::multiprecision::exp(exp_arg) - 1.0;

#if defined(TRAPEZOIDAL_RULE)
  // trapezoidal rule
  const Real x_min = 0.0, &x_max = kTwoPi;
  const Real dx = (x_max - x_min) / number_of_grid_points_;

  normalization_constant_ = 0.0;
  for (int i = 0; i < number_of_grid_points_; ++i)
  {
    normalization_constant_ += (exponent_at_decomposed_point_[i]
        * (integral_of_inv_exp_02pi + periodicity_constraint * integral_up_to_grid_point_[i])
        + exponent_at_decomposed_point_[i + 1]
            * (integral_of_inv_exp_02pi + periodicity_constraint * integral_up_to_grid_point_[i + 1])) * (0.5 * dx);
  } // i
  normalization_constant_ = 1.0 / normalization_constant_;
#elif defined(GAUSS_LEGENDRE_QUADRATURE)
  // Gauss-Legendre quadrature
  normalization_constant_ = 0.0;
  for (int i = 0; i < number_of_grid_points_; ++i)
  {
    for (int m = 0; m < number_of_quadrature_points_; ++m)
    {
      normalization_constant_ += abstract_quadrature_weights_[m]
          * (exponent_at_decomposed_point_[m + i * number_of_quadrature_points_] * (integral_of_inv_exp_02pi
              + periodicity_constraint * integral_up_to_decomposed_point_[m + i * number_of_quadrature_points_]));
    } // m
  } // i
  normalization_constant_ = 1.0 / normalization_constant_;
#endif
}

void SkewedUnimodalCircularDensity::IntegrateForRealPart(MultiprecisionReal &real_part) const
{
  // auxiliary variables
  const MultiprecisionReal &integral_of_inv_exp_02pi
      = integral_up_to_grid_point_[integral_up_to_grid_point_.size() - 1];
  static MultiprecisionReal periodicity_constraint = 0.0, exp_arg = 0.0;
  exp_arg = kTwoPi * group_velocity_ * inverse_diffusion_;
  periodicity_constraint = boost::multiprecision::exp(exp_arg) - 1.0;

#if defined(TRAPEZOIDAL_RULE)
  // expectation E(cos(x)) via a trapezoidal rule
  const Real x_min = 0.0, &x_max = kTwoPi;
  const Real dx = (x_max - x_min) / number_of_grid_points_;
//  static Real x_cur = 0.0, x_nxt = 0.0;

  real_part = 0.0;
  for (int i = 0; i < number_of_grid_points_; ++i)
  {
//    x_cur = i * dx;
//    x_nxt = (i + 1) * dx;
    real_part += (exponent_at_decomposed_point_[i]
        * (integral_of_inv_exp_02pi + periodicity_constraint * integral_up_to_grid_point_[i])
        * cos_of_decomposed_point_[i] + exponent_at_decomposed_point_[i + 1]
        * (integral_of_inv_exp_02pi + periodicity_constraint * integral_up_to_grid_point_[i + 1])
        * cos_of_decomposed_point_[i + 1]) * normalization_constant_ * (0.5 * dx);
  } // i
#elif defined(GAUSS_LEGENDRE_QUADRATURE)
  // expectation E(cos(x)) via Gauss-Legendre quadrature
  real_part = 0.0;
  for (int i = 0; i < number_of_grid_points_; ++i)
  {
    for (int m = 0; m < number_of_quadrature_points_; ++m)
    {
      real_part += abstract_quadrature_weights_[m]
          * (exponent_at_decomposed_point_[m + i * number_of_quadrature_points_] * (integral_of_inv_exp_02pi
              + periodicity_constraint * integral_up_to_decomposed_point_[m + i * number_of_quadrature_points_])
              * normalization_constant_ * cos_of_decomposed_point_[m + i * number_of_quadrature_points_]);
    } // m
  } // i
#endif
}

void SkewedUnimodalCircularDensity::IntegrateForImaginaryPart(MultiprecisionReal &imaginary_part) const
{
  // auxiliary variables
  const MultiprecisionReal &integral_of_inv_exp_02pi
      = integral_up_to_grid_point_[integral_up_to_grid_point_.size() - 1];
  static MultiprecisionReal periodicity_constraint = 0.0, exp_arg = 0.0;
  exp_arg = kTwoPi * group_velocity_ * inverse_diffusion_;
  periodicity_constraint = boost::multiprecision::exp(exp_arg) - 1.0;

#if defined(TRAPEZOIDAL_RULE)
  // expectation E(sin(x)) via a trapezoidal rule
  const Real x_min = 0.0, &x_max = kTwoPi;
  const Real dx = (x_max - x_min) / number_of_grid_points_;
//  static Real x_cur = 0.0, x_nxt = 0.0;

  imaginary_part = 0.0;
  for (int i = 0; i < number_of_grid_points_; ++i)
  {
//    x_cur = i * dx;
//    x_nxt = (i + 1) * dx;
    imaginary_part += (exponent_at_decomposed_point_[i]
        * (integral_of_inv_exp_02pi + periodicity_constraint * integral_up_to_grid_point_[i])
        * sin_of_decomposed_point_[i] + exponent_at_decomposed_point_[i + 1]
        * (integral_of_inv_exp_02pi + periodicity_constraint * integral_up_to_grid_point_[i + 1])
        * sin_of_decomposed_point_[i + 1]) * normalization_constant_ * (0.5 * dx);
  } // i
#elif defined(GAUSS_LEGENDRE_QUADRATURE)
  // expectation E(sin(x)) via Gauss-Legendre quadrature
  imaginary_part = 0.0;
  for (int i = 0; i < number_of_grid_points_; ++i)
  {
    for (int m = 0; m < number_of_quadrature_points_; ++m)
    {
      imaginary_part += abstract_quadrature_weights_[m]
          * (exponent_at_decomposed_point_[m + i * number_of_quadrature_points_] * (integral_of_inv_exp_02pi
              + periodicity_constraint * integral_up_to_decomposed_point_[m + i * number_of_quadrature_points_])
              * normalization_constant_ * sin_of_decomposed_point_[m + i * number_of_quadrature_points_]);
    } // m
  } // i
#endif
}