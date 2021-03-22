//
// Created by Nikita Kruk on 2019-06-03.
//

#ifndef FVMSELFCONSISTENTEQUATIONS_DEFINITIONS_HPP
#define FVMSELFCONSISTENTEQUATIONS_DEFINITIONS_HPP

//#define BCS_CLUSTER
//#define MPI_FOR_PARAMETER_SPAN

#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/math/constants/constants.hpp>

typedef double Real;
typedef std::complex<double> Complex;
//typedef boost::multiprecision::number<boost::multiprecision::cpp_dec_float<500>> MultiprecisionReal;
typedef boost::multiprecision::number<boost::multiprecision::mpfr_float_backend<0>> MultiprecisionReal;

const Complex kI = Complex(0.0, 1.0);
const Real kTwoPi = boost::math::constants::two_pi<Real>();

void GenerateGaussLegendrePoints(Real lower_limit, Real upper_limit, int n, std::vector<Real> &x, std::vector<Real> &w);

#endif //FVMSELFCONSISTENTEQUATIONS_DEFINITIONS_HPP
