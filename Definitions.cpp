//
// Created by Nikita Kruk on 01.01.20.
//

#include "Definitions.hpp"

void GenerateGaussLegendrePoints(Real lower_limit, Real upper_limit, int n, std::vector<Real> &x, std::vector<Real> &w)
{
  const double epsilon = 3.0e-11;
  int n_half = 0, j = 0, i = 0;
  Real middle = 0.0, half_length = 0.0, root = 0.0, root_prev = 0.0;
  Real p_1 = 0.0, p_2 = 0.0, p_3 = 0.0, polynomial_derivative = 0.0;

  n_half = (n + 1) / 2;
  middle = 0.5 * (lower_limit + upper_limit);
  half_length = 0.5 * (upper_limit - lower_limit);
  for (i = 0; i < n_half; ++i)
  {
    root = std::cos(M_PI * (i + 0.75) / (n + 0.5));
    do
    {
      p_1 = 1.0;
      p_2 = 0.0;
      for (j = 0; j < n; ++j)
      {
        p_3 = p_2;
        p_2 = p_1;
        // $(j + 1) P_{j + 1} = (2j + 1)xP_j - jP_{j-1}$
        p_1 = ((2.0 * j + 1.0) * root * p_2 - j * p_3) / (j + 1);
      } // j
      polynomial_derivative = n * (root * p_1 - p_2) / (root * root - 1.0);
      root_prev = root;
      root = root_prev - p_1 / polynomial_derivative;
    } while (std::fabs(root - root_prev) > epsilon);
    x[i] = middle - half_length * root;
    x[n - 1 - i] = middle + half_length * root;
    w[i] = 2.0 * half_length / ((1.0 - root * root) * polynomial_derivative * polynomial_derivative);
    w[n - 1 - i] = w[i];
  } // i
}