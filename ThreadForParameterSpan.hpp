//
// Created by Nikita Kruk on 2019-07-04.
//

#ifndef FVMSELFCONSISTENTEQUATIONS_THREADFORPARAMETERSPAN_HPP
#define FVMSELFCONSISTENTEQUATIONS_THREADFORPARAMETERSPAN_HPP

#include "Thread.hpp"

class ThreadForParameterSpan : public Thread
{
 public:

  ThreadForParameterSpan(int argc, char **argv);
  ~ThreadForParameterSpan();

 private:

};

#endif //FVMSELFCONSISTENTEQUATIONS_THREADFORPARAMETERSPAN_HPP
