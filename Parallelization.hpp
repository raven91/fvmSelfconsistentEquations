//
// Created by Nikita Kruk on 2019-07-04.
//

#ifndef FVMSELFCONSISTENTEQUATIONS_PARALLELIZATION_HPP
#define FVMSELFCONSISTENTEQUATIONS_PARALLELIZATION_HPP

#include "Definitions.hpp"

void LaunchParallelSession(int argc, char **argv);
void FinalizeParallelSession();

#endif //FVMSELFCONSISTENTEQUATIONS_PARALLELIZATION_HPP
