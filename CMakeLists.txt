cmake_minimum_required(VERSION 3.14)
project(fvmSelfconsistentEquations)

set(CMAKE_CXX_STANDARD 17)

include_directories(/Users/nikita/Documents/mpich-install/include/ /usr/local/include/)
link_directories(/Users/nikita/Documents/mpich-install/lib/ /usr/local/lib/)
set(CMAKE_OSX_SYSROOT /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk)

set(GCC_COVERAGE_COMPILE_FLAGS "-march=native -O3 -flto")
set(GCC_COVERAGE_LINK_FLAGS "-lmpc -lmpfr -lgmp")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}")

add_executable(fvmSelfconsistentEquations main.cpp
        SolverHomogeneousZeroLag.cpp SolverHomogeneousZeroLag.hpp
        Definitions.hpp
        HashTuple.hpp
        SolverHomogeneousNonzeroLag.cpp SolverHomogeneousNonzeroLag.hpp
        Parallelization.hpp Parallelization.cpp
        Thread.cpp Thread.hpp
        ThreadForParameterSpan.cpp ThreadForParameterSpan.hpp
        SkewedUnimodalCircularDensity.cpp SkewedUnimodalCircularDensity.hpp
        Definitions.cpp)