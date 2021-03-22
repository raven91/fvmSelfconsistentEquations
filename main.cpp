#include "SolverHomogeneousZeroLag.hpp"
#include "SolverHomogeneousNonzeroLag.hpp"
#include "Parallelization.hpp"
#include "Thread.hpp"
#include "ThreadForParameterSpan.hpp"

int main(int argc, char **argv)
{
  LaunchParallelSession(argc, argv);
  {
#if defined(MPI_FOR_PARAMETER_SPAN)
    ThreadForParameterSpan thread(argc, argv);
#else
    Thread thread(argc, argv);
#endif

    MultiprecisionReal::default_precision(1500); // 1500 used for PRL paper

//  SolverHomogeneousZeroLag solver;
//  solver.FindOrderParameter();
    SolverHomogeneousNonzeroLag solver(&thread);
    solver.FindOrderParameter_coupled();
//    solver.FindModeOfDensityFunction();
//    solver.ComputeSkewnessOfDensityFunction();
  }
  FinalizeParallelSession();

  return 0;
}