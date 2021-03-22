//
// Created by Nikita Kruk on 2019-06-03.
//

#include "SolverHomogeneousNonzeroLag.hpp"
#include "SkewedUnimodalCircularDensity.hpp"

#include <cmath>
#include <iostream>
#include <fstream>
#include <cassert>
#include <limits>
#include <algorithm> // std::fill, std::minmax_element
#include <iterator> // std::distance
#include <string>
#include <chrono>

#include <eigen3/Eigen/Dense>

SolverHomogeneousNonzeroLag::SolverHomogeneousNonzeroLag(Thread *thread) :
    thread_(thread),
    number_of_phase_grid_points_(1000),
    n_sigma_(10),
    n_D_phi_(100),
    n_alpha_(10),
    order_parameters_(),
    velocities_(),
    exponent_at_current_point_(number_of_phase_grid_points_ + 1, MultiprecisionReal(0.0)),
    inverse_exponent_at_current_point_(number_of_phase_grid_points_ + 1, MultiprecisionReal(0.0)),
    integral_up_to_current_point_(number_of_phase_grid_points_ + 1, MultiprecisionReal(0.0)),
    mersenne_twister_generator_(std::random_device{}()),
    uniform_real_distribution_(0.0, 1.0)
{
  sigma_min_ = 0.0;
  sigma_max_ = 10.0;
  dsigma_ = (sigma_max_ - sigma_min_) / n_sigma_;
#if defined(MPI_FOR_PARAMETER_SPAN)
  D_phi_min_ = 0.0;//thread_->GetRank() * 0.3;
  D_phi_max_ = 0.5;//(thread_->GetRank() + 1) * 0.3;
  dD_phi_ = 0.0025;//(D_phi_max_ - D_phi_min_) / n_D_phi_;
  n_D_phi_ = std::nearbyint((D_phi_max_ - D_phi_min_) / dD_phi_);
  dalpha_ = 0.01;
  n_alpha_ = 8;//std::nearbyint((alpha_max_ - alpha_min_) / dalpha_);
  alpha_min_ = 0.0 + thread_->GetRank() * n_alpha_ * dalpha_;//0.0;
  alpha_max_ = 0.0 + (thread_->GetRank() + 1) * n_alpha_ * dalpha_;//1.57;
#else
  D_phi_min_ = 0.03;
  D_phi_max_ = 0.04;
  dD_phi_ = (D_phi_max_ - D_phi_min_) / n_D_phi_;
  alpha_min_ = 1.2;
  alpha_max_ = 1.5;
  dalpha_ = (alpha_max_ - alpha_min_) / n_alpha_;
#endif
}

SolverHomogeneousNonzeroLag::~SolverHomogeneousNonzeroLag()
{

}

//void SolverHomogeneousNonzeroLag::FindOrderParameter_partitive()
//{
//#if defined(__linux__) && defined(BCS_CLUSTER)
//  std::string folder("/home/nkruk/cpp/fvmSelfconsistentEquations/output/HomogeneousSolutionsNonzeroLag/");
//#elif defined(__APPLE__)
//  std::string folder("/Users/nikita/Documents/Projects/fvm/fvmSelfconsistentEquations/HomogeneousSolutionsNonzeroLag/");
//#endif
//  std::ofstream
//      order_parameter_file(folder + std::string("order_parameter_magnitude.txt"), std::ios::out | std::ios::trunc);
//  assert(order_parameter_file.is_open());
//  std::ofstream velocity_file(folder + std::string("velocity.txt"), std::ios::out | std::ios::trunc);
//  assert(velocity_file.is_open());
//
////  for (int i_alpha = 0; i_alpha < n_alpha_; ++i_alpha)
//  {
////    Real alpha = alpha_min_ + i_alpha * dalpha_;
//    Real alpha = 1.5;
//    for (int i_sigma = 0; i_sigma < n_sigma_; ++i_sigma)
//    {
//      Real sigma = sigma_min_ + i_sigma * dsigma_;
////      Real sigma = 6.2;
//      for (int i_D_phi = 0; i_D_phi < n_D_phi_; ++i_D_phi)
//      {
//        Real D_phi = D_phi_min_ + i_D_phi * dD_phi_;
////        Real D_phi = 0.05;
//
//        Real v = -sigma * std::sin(alpha);
//        Real R_prev = 0.5, R_next = 0.5; // initial_guesses
//        bool velocity_successfully_found = false, velocity_is_switched = false;
//        Real velocity_update_fail_counter = 0.0;
//        do // find order parameter
//        {
//          R_prev = R_next;
//
//          /*Real v_min = -sigma * std::sin(alpha), v_max = R_prev;
//          while (std::fabs(v_max - v_min) > 1e-04) // find velocity
//          {
//            v = 0.5 * (v_min + v_max);
//            CalculateAuxiliaryExponentsAndIntegrals(v, R_prev, D_phi, sigma, alpha);
//            Real normalization = ComputeDensityFunctionNormalizationConstant(v, R_prev, D_phi, sigma, alpha);
////            Real normalization_min = ComputeDensityFunctionNormalizationConstant(v, R_prev, D_phi, sigma, alpha);
////            Real normalization_max = ComputeDensityFunctionNormalizationConstant(v, R_prev, D_phi, sigma, alpha);
//            // generally, normalization term is calculated for each velocity separately,
//            // but we need only the sign of the integrals
//            // thus, we assume normalization approx equal for all three cases
//            Real integral = ImaginaryPart(v, R_prev, D_phi, sigma, alpha, normalization);
//            CalculateAuxiliaryExponentsAndIntegrals(v_min, R_prev, D_phi, sigma, alpha);
//            Real integral_min = ImaginaryPart(v_min, R_prev, D_phi, sigma, alpha, normalization);
//            CalculateAuxiliaryExponentsAndIntegrals(v_max, R_prev, D_phi, sigma, alpha);
//            Real integral_max = ImaginaryPart(v_max, R_prev, D_phi, sigma, alpha, normalization);
//            if (integral_min * integral_max < 0.0)
//            {
//              if (integral * integral_min < 0.0)
//              {
//                v_max = v;
//              } else
//              {
//                v_min = v;
//              }
//            } else // if the solution is not in the interval
//            {
//              // TODO: probable infinite cycle
//              break;
////              if (std::fabs(integral_min) < std::fabs(integral_max))
////              {
////                v_max = v;
////                v_min -= std::fabs(v_min);
////              } else
////              {
////                v_min = v;
////                v_max += std::fabs(v_max);
////              }
////              R_prev = uniform_real_distribution_(mersenne_twister_generator_);
////              break;
//            }
//          } // find velocity*/
//
//          /*Real v_prev = -sigma * std::sin(alpha), v_next = -sigma * std::sin(alpha);
//          do // relaxation method
//          {
//            v_prev = v_next;
//            CalculateAuxiliaryExponentsAndIntegrals(v_prev, R_prev, D_phi, sigma, alpha);
//            Real normalization = ComputeDensityFunctionNormalizationConstant(v_prev, R_prev, D_phi, sigma, alpha);
//            Real integral = ImaginaryPart(v_prev, R_prev, D_phi, sigma, alpha, normalization);
//
//            v_next = v_prev + 1e-04;
//            CalculateAuxiliaryExponentsAndIntegrals(v_next, R_prev, D_phi, sigma, alpha);
//            normalization = ComputeDensityFunctionNormalizationConstant(v_next, R_prev, D_phi, sigma, alpha);
//            Real integral_next = ImaginaryPart(v_next, R_prev, D_phi, sigma, alpha, normalization);
//            Real derivative = (integral_next - integral) / 1e-04;
//            if (derivative < 0.0)
//            {
//              v_next = v_prev + 1e-4 * integral;
//            } else
//            {
//              v_next = v_prev - 1e-4 * integral;
//            }
//          } while (std::fabs(v_next - v_prev) > 1e-04);
//          v = v_next;*/
//
//          Real v_next = 0.0;
////          if (!velocity_is_switched)
////          {
//          v_next = -sigma * std::sin(alpha);
//          velocity_successfully_found = false;
//          velocity_update_fail_counter = 0.0;
//          velocity_is_switched = false;
////          }
//          do // Newton's method
//          {
//            v = v_next;
//
//            Real v_left = v - 1e-06;
//            CalculateAuxiliaryExponentsAndIntegrals(v_left, R_prev, D_phi, sigma, alpha);
//            Real normalization = ComputeDensityFunctionNormalizationConstant(v_left, R_prev, D_phi, sigma, alpha);
//            Real integral_left = ImaginaryPart(v_left, R_prev, D_phi, sigma, alpha, normalization);
//
//            CalculateAuxiliaryExponentsAndIntegrals(v, R_prev, D_phi, sigma, alpha);
//            normalization = ComputeDensityFunctionNormalizationConstant(v, R_prev, D_phi, sigma, alpha);
//            Real integral = ImaginaryPart(v, R_prev, D_phi, sigma, alpha, normalization);
//
//            Real derivative = (integral - integral_left) / (v - v_left);
//            Real v_trial = v
//                - std::pow(10.0, -velocity_update_fail_counter) * integral / derivative;
////            Real v_trial = v - (v - v_left) / (integral - integral_left) * integral;
//            if (!velocity_is_switched && (v_trial < -sigma * std::sin(alpha)))
//            {
//              velocity_is_switched = true;
//              velocity_successfully_found = false;
//              velocity_update_fail_counter = 0.0;
//              v_next = 0.0;
//              continue;
//            }
//
//            if (!std::isfinite(v_trial))
//            {
//              break;
//            }
//
//            // velocity can only be negative
//            v_next = std::min(0.0, v_trial);
//            if (velocity_is_switched)
//            {
//              if (derivative < 0.0)
//              {
//                velocity_successfully_found = false;
////                break;
//                v_next = 0.0;
//                velocity_update_fail_counter += 1.0;
//              } else if (std::fabs(v_next) > 1e-04)
//              {
//                velocity_successfully_found = false;
//                velocity_update_fail_counter += 1.0;
//              }
//            }
////            if (velocity_is_switched && (derivative < 0.0))
////            {
////              velocity_successfully_found = false;
////              break;
////            }
////            if (velocity_is_switched && (std::fabs(v_next) > 1e-04))
////            {
////              velocity_successfully_found = false;
////              velocity_update_fail_counter += 1.0;
//////              break;
////            }
//          } while (std::fabs(v_next - v) > 1e-04);
//          if (std::fabs(v_next - v) <= 1e-04)
//          {
//            velocity_successfully_found = true;
//            velocity_update_fail_counter = 0.0;
//            velocity_is_switched = false;
//            v = v_next;
//          }
//
//          CalculateAuxiliaryExponentsAndIntegrals(v, R_prev, D_phi, sigma, alpha);
//          Real normalization = ComputeDensityFunctionNormalizationConstant(v, R_prev, D_phi, sigma, alpha);
//          R_next = std::max(0.0,
//                            RealPart(v,
//                                     R_prev,
//                                     D_phi,
//                                     sigma,
//                                     alpha,
//                                     normalization)); // order parameter can only be in [0,1]
//
//          if (!(std::isfinite(v) && std::isfinite(R_next)))
//          {
//            v = std::numeric_limits<Real>::quiet_NaN();
//            R_next = std::numeric_limits<Real>::quiet_NaN();
//            break;
//          }
//        } while ((std::fabs(R_next - R_prev) > 1e-04) || (!velocity_successfully_found)); // find order parameter
//
////        order_parameters_[std::make_tuple(sigma, D_phi, alpha)] = R_next;
////        velocities_[std::make_tuple(sigma, D_phi, alpha)] = v;
//        std::cout << "sigma:" << sigma << ", D_phi:" << D_phi << ", alpha:" << alpha << ", R:" << R_next << ", v:" << v
//                  << std::endl;
//        order_parameter_file << sigma << '\t' << D_phi << '\t' << alpha << '\t' << R_next << std::endl;
//        velocity_file << sigma << '\t' << D_phi << '\t' << alpha << '\t' << v << std::endl;
//      } // i_alpha
//    } // i_D_phi
//  } // i_sigma
//
//  order_parameter_file.close();
//  velocity_file.close();
//}

void SolverHomogeneousNonzeroLag::FindOrderParameter_coupled()
{
  std::ofstream order_parameter_file, velocity_file;
  CreateOutputFiles(order_parameter_file, velocity_file);

  const Real tolerance = 1e-6; // note: is implicitly dependent on number_of_phase_grid_points_
  const int max_iteration_counter = 20;
  const Real R_min = 0.0, R_max = 1.0;

  MultiprecisionReal normalization = 0.0;
  MultiprecisionReal real_part = 0.0, imaginary_part = 0.0;
  MultiprecisionReal first_equation = 0.0, second_equation = 0.0;
  MultiprecisionReal first_equation_prev_op = 0.0, first_equation_prev_v = 0.0;
  MultiprecisionReal second_equation_prev_op = 0.0, second_equation_prev_v = 0.0;
  Eigen::Vector2d equation(0.0, 0.0);
  Eigen::Vector2d parameter_cur(0.0, 0.0), parameter_next(0.0, 0.0), system_solution(0.0, 0.0);
//  Real jacobian_11 = 0.0, jacobian_12 = 0.0, jacobian_21 = 0.0, jacobian_22 = 0.0;
  Eigen::Matrix2d jacobian = Eigen::Matrix2d::Zero();
  SkewedUnimodalCircularDensity circular_density;

  std::chrono::time_point<std::chrono::system_clock> timer = std::chrono::system_clock::now();

  Real sigma(0.0), D_phi(0.0), alpha(0.0);
  /*std::ifstream input_parameter_file;
  std::ofstream output_parameter_file;
  CreateParameterFiles(input_parameter_file, output_parameter_file);*/
//  std::vector<std::array<Real, 3>> comparison_to_fvm_parameters;
//  DefineParametersForComparisonToFvm(comparison_to_fvm_parameters);
//  for (int i_D_phi = 1; i_D_phi <= n_D_phi_; ++i_D_phi)
  {
//    D_phi = D_phi_min_ + i_D_phi * dD_phi_;
    D_phi = 0.0075;
//    for (int i_alpha = 1; i_alpha <= n_alpha_; ++i_alpha)
//  for (int i_alpha = n_alpha_; i_alpha >= 1; --i_alpha)
    {
//      alpha = alpha_min_ + i_alpha * dalpha_;
      alpha = 1.555;
//      for (int i_sigma = 1; i_sigma <= n_sigma_; ++i_sigma)
//      while (input_parameter_file >> sigma >> D_phi >> alpha)
//      for (const std::array<Real, 3> &comparison_parameter_set : comparison_to_fvm_parameters)
//      sigma = 1.0;
        /*if (D_phi > 0.5 * sigma * std::cos(alpha))
        {
          continue;
        } else*/
      {
//        sigma = sigma_min_ + i_sigma * dsigma_;
        sigma = 1.0;
//        D_phi = 0.99 * (0.5 * sigma * std::cos(alpha));
//        sigma = comparison_parameter_set[0];
//        D_phi = comparison_parameter_set[1];
//        alpha = comparison_parameter_set[2];

        const Real v_min = -sigma * std::sin(alpha), v_max = 0.0;

        // initial parameter values
        Real R_cur = 0.5 * R_max, R_next = 0.5 * R_max;
        Real v_cur = 0.5 * v_min, v_next = 0.5 * v_min;
        int iteration_counter = 0;
        do
        {
          R_cur = R_next;
          v_cur = v_next;
          const Real R_minus_delta = R_cur - 1e-6;
          const Real v_minus_delta = v_cur - 1e-6;
          std::cout << std::setprecision(std::numeric_limits<Real>::max_digits10) << R_cur << ", "
                    << std::setprecision(std::numeric_limits<Real>::max_digits10) << v_cur << std::endl;

//          CalculateAuxiliaryExponentsAndIntegrals(v_cur, R_cur, D_phi, sigma, alpha);
//          ComputeDensityFunctionNormalizationConstant(v_cur, D_phi, normalization);
//          RealPart(v_cur, D_phi, normalization, real_part);
          circular_density.ReinitializeParameters(sigma, alpha, D_phi, v_cur, R_cur);
          circular_density.IntegrateForRealPart(real_part);
          first_equation = real_part - R_cur;
//          CalculateAuxiliaryExponentsAndIntegrals(v_cur, R_minus_delta, D_phi, sigma, alpha);
//          ComputeDensityFunctionNormalizationConstant(v_cur, D_phi, normalization);
//          RealPart(v_cur, D_phi, normalization, real_part);
          circular_density.ReinitializeParameters(sigma, alpha, D_phi, v_cur, R_minus_delta);
          circular_density.IntegrateForRealPart(real_part);
          first_equation_prev_op = real_part - R_minus_delta;
//          jacobian_11 = (first_equation - first_equation_prev_op) / (R_cur - R_minus_delta);
          jacobian(0, 0) = Real(first_equation - first_equation_prev_op) / (R_cur - R_minus_delta);

//          CalculateAuxiliaryExponentsAndIntegrals(v_minus_delta, R_cur, D_phi, sigma, alpha);
//          ComputeDensityFunctionNormalizationConstant(v_minus_delta, D_phi, normalization);
//          RealPart(v_minus_delta, D_phi, normalization, real_part);
          circular_density.ReinitializeParameters(sigma, alpha, D_phi, v_minus_delta, R_cur);
          circular_density.IntegrateForRealPart(real_part);
          first_equation_prev_v = real_part - R_cur;
//          jacobian_12 = (first_equation - first_equation_prev_v) / (v_cur - v_minus_delta);
          jacobian(0, 1) = Real(first_equation - first_equation_prev_v) / (v_cur - v_minus_delta);

//          CalculateAuxiliaryExponentsAndIntegrals(v_cur, R_cur, D_phi, sigma, alpha);
//          ComputeDensityFunctionNormalizationConstant(v_cur, D_phi, normalization);
//          ImaginaryPart(v_cur, D_phi, normalization, imaginary_part);
          circular_density.ReinitializeParameters(sigma, alpha, D_phi, v_cur, R_cur);
          circular_density.IntegrateForImaginaryPart(imaginary_part);
          second_equation = imaginary_part;
//          CalculateAuxiliaryExponentsAndIntegrals(v_cur, R_minus_delta, D_phi, sigma, alpha);
//          ComputeDensityFunctionNormalizationConstant(v_cur, D_phi, normalization);
//          ImaginaryPart(v_cur, D_phi, normalization, imaginary_part);
          circular_density.ReinitializeParameters(sigma, alpha, D_phi, v_cur, R_minus_delta);
          circular_density.IntegrateForImaginaryPart(imaginary_part);
          second_equation_prev_op = imaginary_part;
//          jacobian_21 = (second_equation - second_equation_prev_op) / (R_cur - R_minus_delta);
          jacobian(1, 0) = Real(second_equation - second_equation_prev_op) / (R_cur - R_minus_delta);

//          CalculateAuxiliaryExponentsAndIntegrals(v_minus_delta, R_cur, D_phi, sigma, alpha);
//          ComputeDensityFunctionNormalizationConstant(v_minus_delta, D_phi, normalization);
//          ImaginaryPart(v_minus_delta, D_phi, normalization, imaginary_part);
          circular_density.ReinitializeParameters(sigma, alpha, D_phi, v_minus_delta, R_cur);
          circular_density.IntegrateForImaginaryPart(imaginary_part);
          second_equation_prev_v = imaginary_part;
//          jacobian_22 = (second_equation - second_equation_prev_v) / (v_cur - v_minus_delta);
          jacobian(1, 1) = Real(second_equation - second_equation_prev_v) / (v_cur - v_minus_delta);

//          Real determinant = jacobian(0, 0) * jacobian(1, 1) - jacobian(0, 1) * jacobian(1, 0);
//          R_next = R_cur -
//              (jacobian(1, 1) * first_equation - jacobian(0, 1) * second_equation) / determinant;
//          v_next = v_cur -
//              (-jacobian(1, 0) * first_equation + jacobian(0, 0) * second_equation) / determinant;

          equation(0) = Real(first_equation);
          equation(1) = Real(second_equation);
          parameter_cur(0) = R_cur;
          parameter_cur(1) = v_cur;
          Eigen::HouseholderQR<Eigen::Matrix2d> qr_decomposition(jacobian);
          system_solution = qr_decomposition.solve(equation);
          parameter_next = parameter_cur - system_solution;
          R_next = parameter_next(0);
          v_next = parameter_next(1);

          // if the parameters give undefined solutions
          if (!(std::isfinite(R_next) && std::isfinite(v_next)) && (iteration_counter == 0))
          {
            v_next = std::numeric_limits<Real>::quiet_NaN();
            R_next = std::numeric_limits<Real>::quiet_NaN();
            break;
          }

          // order parameter magnitude must be [0,1]
          /*if (R_next < R_min)
          {
            R_next = R_min;
          } else if (R_next > R_max)
          {
            R_next = R_max;
          }
          // velocity must be [-\sigma*\sin(\alpha),0] for \alpha\in[0,\pi/2]
          if (v_next < v_min)
          {
            v_next = v_min;
          } else if (v_next > v_max)
          {
            v_next = v_max;
          }*/

          // if order parameter is zero, velocity is not uniquely determinded by the system
          if (IsInRange(R_next, -tolerance, tolerance) && !IsInRange(v_next, v_min - tolerance, v_max + tolerance))
          {
            std::cout << "the cycle is stopped: near zero order parameter is reached" << std::endl;
            break;
          }

          ++iteration_counter;
          if (iteration_counter > max_iteration_counter)
          {
            if (!(std::isfinite(R_next) && std::isfinite(v_next)))
            {
              std::cout << "the cycle is renewed: R=" << R_next << ", v=" << v_next << std::endl;

              std::uniform_real_distribution<Real> order_parameter_distribution(R_min, R_max);
              R_next = order_parameter_distribution(mersenne_twister_generator_);
              R_cur = R_next + 2.0 * tolerance;

              std::uniform_real_distribution<Real> velocity_distribution(v_min, v_max);
              v_next = velocity_distribution(mersenne_twister_generator_);
              v_cur = v_next + 2.0 * tolerance;

              iteration_counter = 0;
            } else if (!(IsInRange(R_next, R_min - tolerance, R_max + tolerance)
                && IsInRange(v_next, v_min - tolerance, v_max + tolerance)))
            {
              std::cout << "the cycle is renewed: R=" << R_next << ", v=" << v_next << std::endl;

              std::uniform_real_distribution<Real> order_parameter_distribution(R_min, R_max);
              R_next = order_parameter_distribution(mersenne_twister_generator_);
              R_cur = R_next + 2.0 * tolerance;

              std::uniform_real_distribution<Real> velocity_distribution(v_min, v_max);
              v_next = velocity_distribution(mersenne_twister_generator_);
              v_cur = v_next + 2.0 * tolerance;

              iteration_counter = 0;
            } else
            {
              std::cout << "the cycle is stopped: maximum number of iterations is reached" << std::endl;
              break;
            }
          }
        } while ((std::fabs(R_next - R_cur) > tolerance) || (std::fabs(v_next - v_cur) > tolerance)
            || (!IsInRange(R_next, R_min - tolerance, R_max + tolerance))
            || (!IsInRange(v_next, v_min - tolerance, v_max + tolerance))); // L^1 norm
        if (IsInRange(R_next, -tolerance, tolerance))
        {
          R_next = 0.0;
          v_next = 0.0;
        }

        std::cout << "sigma:" << sigma << ", D_phi:" << D_phi << ", alpha:" << alpha
                  << ", R:" << std::setprecision(std::numeric_limits<Real>::max_digits10) << R_next << ", v:"
                  << std::setprecision(std::numeric_limits<Real>::max_digits10) << v_next << std::endl;
        order_parameter_file << sigma << '\t' << D_phi << '\t' << alpha << '\t'
                             << std::setprecision(std::numeric_limits<Real>::max_digits10) << R_next << std::endl;
        velocity_file << sigma << '\t' << D_phi << '\t' << alpha << '\t'
                      << std::setprecision(std::numeric_limits<Real>::max_digits10) << v_next << std::endl;

        std::chrono::duration<Real> elapsed_seconds = std::chrono::system_clock::now() - timer;
        std::cout << "time per cycle: " << elapsed_seconds.count() << "s" << std::endl;
        timer = std::chrono::system_clock::now();

//        output_parameter_file << sigma << '\t' << D_phi << '\t' << alpha << '\t'
//                              << std::setprecision(std::numeric_limits<Real>::digits10) << R_next << '\t'
//                              << std::setprecision(std::numeric_limits<Real>::digits10) << v_next << std::endl;
      } // i_sigma
    } // i_D_phi
  } // i_alpha

  order_parameter_file.close();
  velocity_file.close();
}

void SolverHomogeneousNonzeroLag::CreateParameterFiles(std::ifstream &input_parameter_file,
                                                       std::ofstream &output_parameter_file) const
{
  input_parameter_file.open(std::string(
      "/Users/nikita/Documents/Projects/fvm/fvmSelfconsistentEquations/input_parameter_values.txt"), std::ios::in);
  assert(input_parameter_file.is_open());

  output_parameter_file.open(std::string(
      "/Users/nikita/Documents/Projects/fvm/fvmSelfconsistentEquations/output_parameter_values.txt"),
                             std::ios::out | std::ios::trunc);
  assert(output_parameter_file.is_open());
}

void SolverHomogeneousNonzeroLag::CreateOutputFiles(std::ofstream &order_parameter_file,
                                                    std::ofstream &velocity_file) const
{
#if defined(__linux__) && defined(BCS_CLUSTER)
  std::string folder("/home/nkruk/cpp/fvmSelfconsistentEquations/output/HomogeneousSolutionsNonzeroLag/max_digits10/10000/");
#elif defined(__APPLE__)
  std::string folder("/Users/nikita/Documents/Projects/fvm/fvmSelfconsistentEquations/HomogeneousSolutionsNonzeroLag/");
#endif

#if defined(MPI_FOR_PARAMETER_SPAN)
  order_parameter_file.open(folder + std::string("order_parameter_magnitude_nq") + std::to_string(thread_->GetRank()) + std::string(".txt"), std::ios::out | std::ios::trunc);
  velocity_file.open(folder + std::string("velocity_nq") + std::to_string(thread_->GetRank()) + std::string(".txt"), std::ios::out | std::ios::trunc);
#else
  order_parameter_file.open(folder + std::string("order_parameter_magnitude.txt"), std::ios::out | std::ios::trunc);
  velocity_file.open(folder + std::string("velocity.txt"), std::ios::out | std::ios::trunc);
#endif
  assert(order_parameter_file.is_open());
  assert(velocity_file.is_open());
}

void SolverHomogeneousNonzeroLag::FindModeOfDensityFunction()
{
  std::ifstream order_parameter_file, velocity_file;
  CreateInputFiles(order_parameter_file, velocity_file);
  std::ofstream minmax_density_file;
  CreateOutputFile(minmax_density_file, "MinmaxDensity", "minmax_density");

  Real sigma = 0.0, D_phi = 0.0, alpha = 0.0, order_parameter = 0.0, velocity = 0.0;
  MultiprecisionReal gamma = 0.0;
  MultiprecisionReal c_0 = 0.0, c_e = 0.0, I_2pi = 0.0;
  std::vector<MultiprecisionReal> density(number_of_phase_grid_points_ + 1, 0.0);

  while (order_parameter_file >> sigma >> D_phi >> alpha >> order_parameter)
  {
    velocity_file >> sigma >> D_phi >> alpha >> velocity;
    if (!(std::isfinite(order_parameter) && std::isfinite(velocity)))
    {
      continue;
    }
    gamma = sigma * order_parameter / D_phi;
    std::fill(density.begin(), density.end(), MultiprecisionReal(0.0));

    CalculateAuxiliaryExponentsAndIntegrals(velocity, order_parameter, D_phi, sigma, alpha);
    c_e = boost::multiprecision::exp(MultiprecisionReal(kTwoPi * velocity / D_phi)) - MultiprecisionReal(1.0);
    I_2pi = integral_up_to_current_point_[integral_up_to_current_point_.size() - 1];
    DensityFunctionNormalization(velocity, D_phi, c_0);
    for (int i = 0; i < number_of_phase_grid_points_ + 1; ++i)
    {
      density[i] = exponent_at_current_point_[i] * (I_2pi + c_e * integral_up_to_current_point_[i]) / c_0;
    } // i

    const Real phi_min = 0.0, &phi_max = kTwoPi;
    const Real dphi = (phi_max - phi_min) / number_of_phase_grid_points_;
    auto minmax = std::minmax_element(density.begin(), density.end());
    std::cout << sigma << '\t' << D_phi << '\t' << alpha << '\t' << order_parameter << '\t' << velocity << '\t'
              << phi_min + std::distance(density.begin(), minmax.first) * dphi << '\t'
              << phi_min + std::distance(density.begin(), minmax.second) * dphi << std::endl;
    minmax_density_file << sigma << '\t' << D_phi << '\t' << alpha << '\t'
                        << phi_min + std::distance(density.begin(), minmax.first) * dphi << '\t'
                        << phi_min + std::distance(density.begin(), minmax.second) * dphi << std::endl;
  }

  minmax_density_file.close();
  velocity_file.close();
  order_parameter_file.close();
}

void SolverHomogeneousNonzeroLag::ComputeSkewnessOfDensityFunction()
{
  std::ifstream order_parameter_file, velocity_file;
  CreateInputFiles(order_parameter_file, velocity_file);
  std::ofstream skewness_file;
  CreateOutputFile(skewness_file, "Skewness", "circular_skewness");

  Real sigma = 0.0, D_phi = 0.0, alpha = 0.0, order_parameter = 0.0, velocity = 0.0;
  MultiprecisionReal gamma = 0.0;
  MultiprecisionReal c_0 = 0.0, c_e = 0.0, I_2pi = 0.0;
  std::vector<MultiprecisionReal> density(number_of_phase_grid_points_ + 1, 0.0);
  MultiprecisionReal expectation_real(0.0), expectation_imag(0.0),
      mean_direction(0.0), circular_variance(0.0), circular_skewness(0.0);

  while (order_parameter_file >> sigma >> D_phi >> alpha >> order_parameter)
  {
    velocity_file >> sigma >> D_phi >> alpha >> velocity;
    if (!(std::isfinite(order_parameter) && std::isfinite(velocity)))
    {
      continue;
    }
    gamma = sigma * order_parameter / D_phi;
    std::fill(density.begin(), density.end(), MultiprecisionReal(0.0));

    CalculateAuxiliaryExponentsAndIntegrals(velocity, order_parameter, D_phi, sigma, alpha);
    c_e = boost::multiprecision::exp(MultiprecisionReal(kTwoPi * velocity / D_phi)) - MultiprecisionReal(1.0);
    I_2pi = integral_up_to_current_point_[integral_up_to_current_point_.size() - 1];
    DensityFunctionNormalization(velocity, D_phi, c_0);
    for (int i = 0; i < number_of_phase_grid_points_ + 1; ++i)
    {
      density[i] = exponent_at_current_point_[i] * (I_2pi + c_e * integral_up_to_current_point_[i]) / c_0;
    } // i

    const Real phi_min = 0.0, &phi_max = kTwoPi;
    const Real dphi = (phi_max - phi_min) / number_of_phase_grid_points_;
    Real phi_cur = 0.0, phi_nxt = 0.0;
    Complex z_cur(0.0, 0.0), z_nxt(0.0, 0.0);
    expectation_real = 0.0;
    expectation_imag = 0.0;
    for (int i = 0; i < number_of_phase_grid_points_; ++i)
    { // trapezoidal rule
      phi_cur = phi_min + i * dphi;
      phi_nxt = phi_min + (i + 1) * dphi;
      z_cur = std::exp(kI * phi_cur);
      z_nxt = std::exp(kI * phi_nxt);
      expectation_real += (z_cur.real() * density[i] + z_nxt.real() * density[i + 1]) * (0.5 * dphi);
      expectation_imag += (z_cur.imag() * density[i] + z_nxt.imag() * density[i + 1]) * (0.5 * dphi);
    } // i
    mean_direction = boost::multiprecision::atan2(expectation_imag, expectation_real);

    circular_variance =
        1.0 - boost::multiprecision::sqrt(expectation_real * expectation_real + expectation_imag * expectation_imag);
    /*for (int i = 0; i < number_of_phase_grid_points_; ++i)
    { // trapezoidal rule
      phi_cur = phi_min + i * dphi;
      phi_nxt = phi_min + (i + 1) * dphi;
      variance += ((phi_cur - expectation) * (phi_cur - expectation) * density[i]
          + (phi_nxt - expectation) * (phi_nxt - expectation) * density[i + 1]) * (0.5 * dphi);
    } // i
    standard_deviation = boost::multiprecision::sqrt(variance);*/

    circular_skewness = 0.0;
    for (int i = 0; i < number_of_phase_grid_points_; ++i)
    { // trapezoidal rule
      phi_cur = phi_min + i * dphi;
      phi_nxt = phi_min + (i + 1) * dphi;
//      skewness += (boost::multiprecision::pow((phi_cur - expectation) / standard_deviation, 3.0) * density[i]
//          + boost::multiprecision::pow((phi_nxt - expectation) / standard_deviation, 3.0) * density[i + 1])
//          * (0.5 * dphi);
      circular_skewness += (boost::multiprecision::sin(2.0 * (phi_cur - mean_direction)) * density[i]
          + boost::multiprecision::sin(2.0 * (phi_nxt - mean_direction)) * density[i + 1]) * (0.5 * dphi);
    } // i
    circular_skewness /= boost::multiprecision::pow(circular_variance, 1.5);

    std::cout << sigma << '\t' << D_phi << '\t' << alpha << '\t' << order_parameter << '\t' << velocity << '\t'
              << mean_direction << '\t' << circular_variance << '\t' << circular_skewness << std::endl;
    skewness_file << sigma << '\t' << D_phi << '\t' << alpha << '\t'
                  << mean_direction << '\t' << circular_variance << '\t' << circular_skewness << std::endl;
  }

  skewness_file.close();
  velocity_file.close();
  order_parameter_file.close();
}

void SolverHomogeneousNonzeroLag::CreateInputFiles(std::ifstream &order_parameter_file,
                                                   std::ifstream &velocity_file) const
{
#if defined(__linux__) && defined(BCS_CLUSTER)
  std::string input_folder("/home/nkruk/cpp/fvmSelfconsistentEquations/output/HomogeneousSolutionsNonzeroLag/");
#elif defined(__APPLE__)
  std::string
      input_folder("/Users/nikita/Documents/Projects/fvm/fvmSelfconsistentEquations/HomogeneousSolutionsNonzeroLag/");
#endif

#if defined(MPI_FOR_PARAMETER_SPAN)
  order_parameter_file.open(input_folder + std::string("order_parameter_magnitude_s") + std::to_string(thread_->GetRank()) + std::string(".txt"), std::ios::in);
  velocity_file.open(input_folder + std::string("velocity_s") + std::to_string(thread_->GetRank()) + std::string(".txt"), std::ios::in);
#else
  order_parameter_file.open(input_folder + std::string("order_parameter_magnitude.txt"), std::ios::in);
  velocity_file.open(input_folder + std::string("velocity.txt"), std::ios::in);
#endif
  assert(order_parameter_file.is_open());
  assert(velocity_file.is_open());
}

void SolverHomogeneousNonzeroLag::CreateOutputFile(std::ofstream &file,
                                                   const std::string &subfolder,
                                                   const std::string &base_name) const
{
#if defined(__linux__) && defined(BCS_CLUSTER)
  std::string output_folder("/home/nkruk/cpp/fvmSelfconsistentEquations/output/" + subfolder + "/");
#elif defined(__APPLE__)
  std::string output_folder("/Users/nikita/Documents/Projects/fvm/fvmSelfconsistentEquations/" + subfolder + "/");
#endif

#if defined(MPI_FOR_PARAMETER_SPAN)
  file.open(output_folder + std::string(base_name + "_s") + std::to_string(thread_->GetRank())
                                        + std::string(".txt"), std::ios::out | std::ios::trunc);
#else
  file.open(output_folder + std::string(base_name + ".txt"), std::ios::out | std::ios::trunc);
#endif
  assert(file.is_open());
}

void SolverHomogeneousNonzeroLag::DefineParametersForComparisonToFvm(std::vector<std::array<Real,
                                                                                            3>> &comparison_to_fvm_parameters) const
{
  Real sigma = 1.0, alpha = 0.0, D_phi = 0.0;

  alpha = 0.5;
  for (int i = 0; i < 200; ++i)
  {
    D_phi = 0.0025 + 0.0025 * i;
    comparison_to_fvm_parameters.push_back({sigma, D_phi, alpha});
  } // i
  alpha = 1.0;
  for (int i = 0; i < 132; ++i)
  {
    D_phi = 0.0025 + 0.0025 * i;
    comparison_to_fvm_parameters.push_back({sigma, D_phi, alpha});
  } // i
  alpha = 1.5;
  for (int i = 0; i < 36; ++i)
  {
    D_phi = 0.0025 + 0.0025 * i;
    comparison_to_fvm_parameters.push_back({sigma, D_phi, alpha});
  } // i

  D_phi = 0.1;
  for (int i = 0; i < 155; ++i)
  {
    alpha = 0.0 + 0.01 * i;
    comparison_to_fvm_parameters.push_back({sigma, D_phi, alpha});
  } // i
  D_phi = 0.25;
  for (int i = 0; i < 120; ++i)
  {
    alpha = 0.0 + 0.01 * i;
    comparison_to_fvm_parameters.push_back({sigma, D_phi, alpha});
  } // i
  D_phi = 0.4;
  for (int i = 0; i < 81; ++i)
  {
    alpha = 0.0 + 0.01 * i;
    comparison_to_fvm_parameters.push_back({sigma, D_phi, alpha});
  } // i
}

//// $ E(\varphi,v) = e^{-v / D_\varphi * \varphi + \sigma / D_\varphi * R * \cos(\varphi + \alpha)} $
//Real SolverHomogeneousNonzeroLag::E(Real phi, Real v, Real R, Real D_phi, Real sigma, Real alpha) const
//{
//  return std::exp((-v * phi + sigma * R * std::cos(phi + alpha)) / D_phi);
//}

// $ E(\varphi,v) = e^{-v / D_\varphi * \varphi + \sigma / D_\varphi * R * \cos(\varphi + \alpha)} $
void SolverHomogeneousNonzeroLag::E(Real phi,
                                    Real v,
                                    Real R,
                                    Real D_phi,
                                    Real sigma,
                                    Real alpha,
                                    MultiprecisionReal &res) const
{
  static MultiprecisionReal exp_arg = 0.0;
  exp_arg = (-v * phi + sigma * R * std::cos(phi + alpha)) / D_phi;
  res = boost::multiprecision::exp(exp_arg);
}

//// $ E(\varphi,v) = e^{v / D_\varphi * \varphi - \sigma / D_\varphi * R * \cos(\varphi + \alpha)} $
//Real SolverHomogeneousNonzeroLag::EInv(Real phi, Real v, Real R, Real D_phi, Real sigma, Real alpha) const
//{
//  return std::exp((v * phi - sigma * R * std::cos(phi + alpha)) / D_phi);
//}

// $ E(\varphi,v) = e^{v / D_\varphi * \varphi - \sigma / D_\varphi * R * \cos(\varphi + \alpha)} $
void SolverHomogeneousNonzeroLag::EInv(Real phi,
                                       Real v,
                                       Real R,
                                       Real D_phi,
                                       Real sigma,
                                       Real alpha,
                                       MultiprecisionReal &res) const
{
  static MultiprecisionReal exp_arg = 0.0;
  exp_arg = (v * phi - sigma * R * std::cos(phi + alpha)) / D_phi;
  res = boost::multiprecision::exp(exp_arg);
}

//// $ IntegrateInverseExponentUpTo(\varphi,v) = \int_0^\varphi E^{-1}(\varphi',v) \mathrm{d}\varphi' $
//// $ IntegrateInverseExponentUpTo(\varphi,v) = \int_0^\varphi e^{v / D_\varphi * \varphi - \sigma / D_\varphi * R * \cos(\varphi + \alpha)} \mathrm{d}\varphi' $
//Real SolverHomogeneousNonzeroLag::IntegrateInverseExponentUpTo(Real phi,
//                                                                  Real v,
//                                                                  Real R,
//                                                                  Real D_phi,
//                                                                  Real sigma,
//                                                                  Real alpha) const
//{
//  Real omega_min = 0, omega_max = phi;
//  int n_omega = number_of_phase_grid_points_;
//  Real domega = (omega_max - omega_min) / n_omega;
//  Real integral = 0.0;
//  Real omega_cur = 0.0, omega_nxt = 0.0;
//
//  // trapezoidal rule
//  for (int i = 0; i < n_omega; ++i)
//  {
//    omega_cur = omega_min + i * domega;
//    omega_nxt = omega_min + (i + 1) * domega;
//    integral +=
//        (EInv(omega_cur, v, R, D_phi, sigma, alpha) + EInv(omega_nxt, v, R, D_phi, sigma, alpha)) * 0.5 * domega;
////    integral += (inverse_exponent_at_current_point_[i] + inverse_exponent_at_current_point_[i + 1]) * 0.5 * domega;
//  } // i
//
//  return integral;
//}

//// $ IntegrateInverseExponentUpTo(\varphi,v) = \int_0^{2\pi} E^{-1}(\varphi',v) \mathrm{d}\varphi' $
//// $ IntegrateInverseExponentUpTo(\varphi,v) = \int_0^{2\pi} e^{v / D_\varphi * \varphi - \sigma / D_\varphi * R * \cos(\varphi + \alpha)} \mathrm{d}\varphi' $
//void SolverHomogeneousNonzeroLag::IntegrateInverseExponentOverRange(Real lower_limit,
//                                                                       Real upper_limit,
//                                                                       Real v,
//                                                                       Real R,
//                                                                       Real D_phi,
//                                                                       Real sigma,
//                                                                       Real alpha)
//{
//  int n_omega = number_of_phase_grid_points_;
//  Real domega = (upper_limit - lower_limit) / n_omega;
//
//  integral_up_to_current_point_[0] = 0.0;
//  for (int i = 0; i < n_omega; ++i)
//  {
//    integral_up_to_current_point_[i + 1] = integral_up_to_current_point_[i]
//        + (inverse_exponent_at_current_point_[i] + inverse_exponent_at_current_point_[i + 1]) * 0.5 * domega;
//  } // i
//}

void SolverHomogeneousNonzeroLag::IntegrateInverseExponentOverRange(Real lower_limit, Real upper_limit)
{
  const Real domega = (upper_limit - lower_limit) / number_of_phase_grid_points_;

  integral_up_to_current_point_[0] = 0.0;
  for (int i = 0; i < number_of_phase_grid_points_; ++i)
  {
    integral_up_to_current_point_[i + 1] = integral_up_to_current_point_[i]
        + (inverse_exponent_at_current_point_[i] + inverse_exponent_at_current_point_[i + 1])
            * (0.5 * domega);
  } // i
}

/*Real SolverHomogeneousNonzeroLag::ComputeDensityFunctionNormalizationConstant(Real v,
                                                                  Real R,
                                                                  Real D_phi,
                                                                  Real sigma,
                                                                  Real alpha) const
{
  Real omega_min = 0, omega_max = 2.0 * M_PI;
  int n_omega = number_of_phase_grid_points_;
  Real domega = (omega_max - omega_min) / n_omega;
//  Real omega_cur = 0.0, omega_nxt = 0.0;

  Real integral = 0.0;
  // auxiliary variables
  Real I_2pi = integral_up_to_current_point_[integral_up_to_current_point_.size() - 1];
  Real c_e = std::exp(2.0 * M_PI * v / D_phi);
  Real fraction_of_integrals = 0.0, fraction_of_integrals_next = 0.0;
  // trapezoidal rule
  for (int i = 0; i < n_omega; ++i)
  {
//    omega_cur = omega_min + i * domega;
//    omega_nxt = omega_min + (i + 1) * domega;
//    integral +=
//        (E(omega_cur, v, R, D_phi, sigma, alpha) * std::fabs(I_2pi + c_e * IntegrateInverseExponentUpTo(omega_cur, v, R, D_phi, sigma, alpha))
//            + E(omega_nxt, v, R, D_phi, sigma, alpha)
//                * std::fabs(I_2pi + c_e * IntegrateInverseExponentUpTo(omega_nxt, v, R, D_phi, sigma, alpha))) * 0.5 * domega;
    fraction_of_integrals = integral_up_to_current_point_[i] / I_2pi;
    fraction_of_integrals_next = integral_up_to_current_point_[i + 1] / I_2pi;
    integral += (exponent_at_current_point_[i] * std::fabs((1.0 - fraction_of_integrals)
                                                               + c_e * fraction_of_integrals)
        + exponent_at_current_point_[i + 1] * std::fabs((1.0 - fraction_of_integrals_next)
                                                            + c_e * fraction_of_integrals_next))
        * 0.5 * domega;
//    integral += (std::exp(std::log(exponent_at_current_point_[i]) + std::log(std::fabs(I_2pi + c_e * integral_up_to_current_point_[i])))
//        + std::exp(std::log(exponent_at_current_point_[i + 1]) + std::log(std::fabs(I_2pi + c_e * integral_up_to_current_point_[i + 1])))) * 0.5
//        * domega;
  } // i

  return integral;
}*/

void SolverHomogeneousNonzeroLag::DensityFunctionNormalization(Real v,
                                                               Real D_phi,
                                                               MultiprecisionReal &normalization) const
{
  const Real omega_min = 0.0, &omega_max = kTwoPi;
  const Real domega = (omega_max - omega_min) / number_of_phase_grid_points_;

  normalization = 0.0;
  // auxiliary variables
  const MultiprecisionReal &I_2pi = integral_up_to_current_point_[integral_up_to_current_point_.size() - 1];
  static MultiprecisionReal c_e = 0.0, exp_arg = 0.0;
  exp_arg = kTwoPi * v / D_phi;
  c_e = boost::multiprecision::exp(exp_arg) - 1.0;

  // trapezoidal rule
  for (int i = 0; i < number_of_phase_grid_points_; ++i)
  {
    normalization += (exponent_at_current_point_[i] * (I_2pi + c_e * integral_up_to_current_point_[i])
        + exponent_at_current_point_[i + 1] * (I_2pi + c_e * integral_up_to_current_point_[i + 1]))
        * (0.5 * domega);
  } // i
}

/*Real SolverHomogeneousNonzeroLag::ImaginaryPart(Real v,
                                                   Real R,
                                                   Real D_phi,
                                                   Real sigma,
                                                   Real alpha,
                                                   Real normalization) const
{
  Real phi_min = 0, phi_max = 2.0 * M_PI;
  int n_phi = number_of_phase_grid_points_;
  Real dphi = (phi_max - phi_min) / n_phi;
  Real phi_cur = 0.0, phi_nxt = 0.0;

  Real integral = 0.0;
  // auxiliary variables
  Real I_2pi = integral_up_to_current_point_[integral_up_to_current_point_.size() - 1];
  Real c_e = std::exp(2.0 * M_PI * v / D_phi);
  Real fraction_of_integrals = 0.0, fraction_of_integrals_next = 0.0;
  // trapezoidal rule
  for (int i = 0; i < n_phi; ++i)
  {
    phi_cur = phi_min + i * dphi;
    phi_nxt = phi_min + (i + 1) * dphi;
//    integral += (E(phi_cur, v, R, D_phi, sigma, alpha) * std::fabs(I_2pi + c_e * IntegrateInverseExponentUpTo(phi_cur, v, R, D_phi, sigma, alpha))
//        * std::sin(phi_cur)
//        + E(phi_nxt, v, R, D_phi, sigma, alpha) * std::fabs(I_2pi + c_e * IntegrateInverseExponentUpTo(phi_nxt, v, R, D_phi, sigma, alpha))
//            * std::sin(phi_nxt)) / normalization * 0.5 * dphi;
    fraction_of_integrals = integral_up_to_current_point_[i] / I_2pi;
    fraction_of_integrals_next = integral_up_to_current_point_[i + 1] / I_2pi;
    integral += (exponent_at_current_point_[i] * std::fabs((1.0 - fraction_of_integrals)
                                                               + c_e * fraction_of_integrals) * std::sin(phi_cur)
        + exponent_at_current_point_[i + 1] * std::fabs((1.0 - fraction_of_integrals_next)
                                                            + c_e * fraction_of_integrals_next) * std::sin(phi_nxt))
        / normalization * 0.5 * dphi;
  } // i

  return integral;
}*/

void SolverHomogeneousNonzeroLag::ImaginaryPart(Real v,
                                                Real D_phi,
                                                const MultiprecisionReal &normalization,
                                                MultiprecisionReal &imaginary_part) const
{
  const Real phi_min = 0.0, &phi_max = kTwoPi;
  const Real dphi = (phi_max - phi_min) / number_of_phase_grid_points_;
  Real phi_cur = 0.0, phi_nxt = 0.0;

  imaginary_part = 0.0;
  // auxiliary variables
  const MultiprecisionReal &I_2pi = integral_up_to_current_point_[integral_up_to_current_point_.size() - 1];
  static MultiprecisionReal c_e = 0.0, exp_arg = 0.0;
  exp_arg = kTwoPi * v / D_phi;
  c_e = boost::multiprecision::exp(exp_arg) - 1.0;

  // trapezoidal rule
  for (int i = 0; i < number_of_phase_grid_points_; ++i)
  {
    phi_cur = phi_min + i * dphi;
    phi_nxt = phi_min + (i + 1) * dphi;
    imaginary_part +=
        (exponent_at_current_point_[i] * (I_2pi + c_e * integral_up_to_current_point_[i]) * std::sin(phi_cur)
            + exponent_at_current_point_[i + 1] * (I_2pi + c_e * integral_up_to_current_point_[i + 1])
                * std::sin(phi_nxt)) / normalization * (0.5 * dphi);
  } // i
}

/*Real SolverHomogeneousNonzeroLag::RealPart(Real v,
                                              Real R,
                                              Real D_phi,
                                              Real sigma,
                                              Real alpha,
                                              Real normalization) const
{
  Real phi_min = 0, phi_max = 2.0 * M_PI;
  int n_phi = number_of_phase_grid_points_;
  Real dphi = (phi_max - phi_min) / n_phi;
  Real phi_cur = 0.0, phi_nxt = 0.0;

  Real integral = 0.0;
  // auxiliary variables
  Real I_2pi = integral_up_to_current_point_[integral_up_to_current_point_.size() - 1];
  Real c_e = std::exp(2.0 * M_PI * v / D_phi);
  Real fraction_of_integrals = 0.0, fraction_of_integrals_next = 0.0;
  // trapezoidal rule
  for (int i = 0; i < n_phi; ++i)
  {
    phi_cur = phi_min + i * dphi;
    phi_nxt = phi_min + (i + 1) * dphi;
//    integral += (E(phi_cur, v, R, D_phi, sigma, alpha) * std::fabs(I_2pi + c_e * IntegrateInverseExponentUpTo(phi_cur, v, R, D_phi, sigma, alpha))
//        * std::cos(phi_cur)
//        + E(phi_nxt, v, R, D_phi, sigma, alpha) * std::fabs(I_2pi + c_e * IntegrateInverseExponentUpTo(phi_nxt, v, R, D_phi, sigma, alpha))
//            * std::cos(phi_nxt)) / normalization * 0.5 * dphi;
    fraction_of_integrals = integral_up_to_current_point_[i] / I_2pi;
    fraction_of_integrals_next = integral_up_to_current_point_[i + 1] / I_2pi;
    integral += (exponent_at_current_point_[i] * std::fabs((1.0 - fraction_of_integrals)
                                                               + c_e * fraction_of_integrals) * std::cos(phi_cur)
        + exponent_at_current_point_[i + 1] * std::fabs((1.0 - fraction_of_integrals_next)
                                                            + c_e * fraction_of_integrals_next) * std::cos(phi_nxt))
        / normalization * 0.5 * dphi;
  } // i

  return integral;
}*/

void SolverHomogeneousNonzeroLag::RealPart(Real v,
                                           Real D_phi,
                                           const MultiprecisionReal &normalization,
                                           MultiprecisionReal &real_part) const
{
  const Real phi_min = 0.0, &phi_max = kTwoPi;
  const Real dphi = (phi_max - phi_min) / number_of_phase_grid_points_;
  Real phi_cur = 0.0, phi_nxt = 0.0;

  real_part = 0.0;
  // auxiliary variables
  const MultiprecisionReal &I_2pi = integral_up_to_current_point_[integral_up_to_current_point_.size() - 1];
  static MultiprecisionReal c_e = 0.0, exp_arg = 0.0;
  exp_arg = kTwoPi * v / D_phi;
  c_e = boost::multiprecision::exp(exp_arg) - 1.0;

  // trapezoidal rule
  for (int i = 0; i < number_of_phase_grid_points_; ++i)
  {
    phi_cur = phi_min + i * dphi;
    phi_nxt = phi_min + (i + 1) * dphi;
    real_part += (exponent_at_current_point_[i] * (I_2pi + c_e * integral_up_to_current_point_[i]) * std::cos(phi_cur)
        + exponent_at_current_point_[i + 1] * (I_2pi + c_e * integral_up_to_current_point_[i + 1]) * std::cos(phi_nxt))
        / normalization * (0.5 * dphi);
  } // i
}

void SolverHomogeneousNonzeroLag::CalculateAuxiliaryExponentsAndIntegrals(Real v,
                                                                          Real R,
                                                                          Real D_phi,
                                                                          Real sigma,
                                                                          Real alpha)
{
  const Real &upper_limit = kTwoPi, lower_limit = 0.0;
  const Real domega = (upper_limit - lower_limit) / number_of_phase_grid_points_;
  Real omega_i = 0.0;
  for (int i = 0; i < number_of_phase_grid_points_ + 1; ++i)
  {
    omega_i = i * domega;
    E(omega_i, v, R, D_phi, sigma, alpha, exponent_at_current_point_[i]);
    EInv(omega_i, v, R, D_phi, sigma, alpha, inverse_exponent_at_current_point_[i]);
  } // i

  IntegrateInverseExponentOverRange(lower_limit, upper_limit);
}

bool SolverHomogeneousNonzeroLag::IsInRange(Real p, Real p_min, Real p_max) const
{
  return (p >= p_min && p <= p_max);
}