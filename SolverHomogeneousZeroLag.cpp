//
// Created by Nikita Kruk on 2019-06-03.
//

#include "SolverHomogeneousZeroLag.hpp"

#include <cmath>
#include <boost/math/special_functions/bessel.hpp>
#include <iostream>
#include <fstream>
#include <cassert>
#include <limits>

SolverHomogeneousZeroLag::SolverHomogeneousZeroLag(int n_sigma, int n_D_phi) :
    n_sigma_(n_sigma),
    n_D_phi_(n_D_phi),
    order_parameters_()
{
  sigma_min_ = 0.0;
  sigma_max_ = 5.0;
  dsigma_ = (sigma_max_ - sigma_min_) / n_sigma_;
  D_phi_min_ = 0.0;
  D_phi_max_ = 0.5;
  dD_phi_ = (D_phi_max_ - D_phi_min_) / n_D_phi_;
}

SolverHomogeneousZeroLag::~SolverHomogeneousZeroLag()
{

}

void SolverHomogeneousZeroLag::FindOrderParameter()
{
  const Real tolerance = 1e-15;
//  for (int i_sigma = 1; i_sigma <= n_sigma_; ++i_sigma)
  {
//    const Real sigma = sigma_min_ + i_sigma * dsigma_;
    const Real sigma = 1.0;
    for (int i_D_phi = 1; i_D_phi <= n_D_phi_; ++i_D_phi)
    {
      const Real D_phi = D_phi_min_ + i_D_phi * dD_phi_;
//      const Real D_phi = 0.43;
      Real R_prev = 0.5, R_next = 0.5; // initial guesses
      do
      {
        R_prev = R_next;
        try
        {
          R_next = Real(boost::math::cyl_bessel_i(1, sigma * R_prev / (D_phi))
                            / boost::math::cyl_bessel_i(0, sigma * R_prev / (D_phi)));
        } catch (const std::exception &e)
        {
          break;
        }
        if (!std::isfinite(R_next))
        {
//          R_next = std::numeric_limits<Real>::quiet_NaN();
          break;
        }
      } while (std::fabs(R_next - R_prev) > tolerance);
      order_parameters_[std::make_tuple(sigma, D_phi)] = R_next;
      std::cout << "sigma:" << sigma << ", D_phi:" << D_phi << ", R:"
                << std::setprecision(std::numeric_limits<Real>::max_digits10) << R_next << std::endl;
    } // i_D_phi
  } // i_sigma

#if defined(__linux__) && defined(BCS_CLUSTER)
  std::string folder("/home/nkruk/cpp/fvmSelfconsistentEquations/output/HomogeneousSolutionsZeroLag/");
#elif defined(__APPLE__)
  std::string folder("/Users/nikita/Documents/Projects/fvm/fvmSelfconsistentEquations/HomogeneousSolutionsZeroLag/");
#endif
  std::ofstream
      order_parameter_file(folder + std::string("order_parameter_magnitude.txt"), std::ios::out | std::ios::trunc);
  assert(order_parameter_file.is_open());
  for (auto op : order_parameters_)
  {
    order_parameter_file << std::get<0>(op.first) << '\t' << std::get<1>(op.first) << '\t'
                         << std::setprecision(std::numeric_limits<Real>::max_digits10) << op.second << std::endl;
  } // op
  order_parameter_file.close();
}