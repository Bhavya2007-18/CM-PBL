cmake_minimum_required(VERSION 3.20)
project(EngineeringEquationSolver VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(MSVC)
    add_compile_options(/W4 /EHsc /utf-8)
else()
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# Core Math & Infrastructure Library
set(EES_CORE_SOURCES
    src/core/Vector.cpp
    src/core/Matrix.cpp
    src/core/ExpressionParser.cpp
    src/core/ErrorAnalysis.cpp
    src/core/NumericalResult.cpp
)
add_library(ees_core STATIC ${EES_CORE_SOURCES})
target_include_directories(ees_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Unit I Library (Foundations, Roots, Optimization)
set(EES_UNIT1_SOURCES
    src/unit1/Foundations.cpp
    src/unit1/RootFinding.cpp
    src/unit1/Optimization.cpp
)
add_library(ees_unit1 STATIC ${EES_UNIT1_SOURCES})
target_link_libraries(ees_unit1 PUBLIC ees_core)

# Unit II Library (Interpolation, Integration)
set(EES_UNIT2_SOURCES
    src/unit2/Interpolation.cpp
    src/unit2/Integration.cpp
)
add_library(ees_unit2 STATIC ${EES_UNIT2_SOURCES})
target_link_libraries(ees_unit2 PUBLIC ees_core)

# Unit III Library (Linear Systems, Factorization, Eigenvalues, Splines)
set(EES_UNIT3_SOURCES
    src/unit3/LinearSystems.cpp
    src/unit3/Splines.cpp
)
add_library(ees_unit3 STATIC ${EES_UNIT3_SOURCES})
target_link_libraries(ees_unit3 PUBLIC ees_core)

# Unit IV Library (ODEs, PDEs)
set(EES_UNIT4_SOURCES
    src/unit4/ODESolvers.cpp
    src/unit4/PDESolvers.cpp
)
add_library(ees_unit4 STATIC ${EES_UNIT4_SOURCES})
target_link_libraries(ees_unit4 PUBLIC ees_core)

# Application Libraries (Comparison, Problems, Advisor, UI, CLI)
set(EES_APP_SOURCES
    src/comparison/MethodComparator.cpp
    src/problems/ProblemLibrary.cpp
    src/advisor/MethodAdvisor.cpp
    src/ui/Visualizer.cpp
    src/ui/Workbench.cpp
    src/ui/CLIParser.cpp
)
add_library(ees_workbench STATIC ${EES_APP_SOURCES})
target_link_libraries(ees_workbench PUBLIC ees_core ees_unit1 ees_unit2 ees_unit3 ees_unit4)

# Main Executable Targets: ees and ees_app
add_executable(ees src/main.cpp)
target_link_libraries(ees PRIVATE ees_workbench)

add_executable(ees_app src/main.cpp)
target_link_libraries(ees_app PRIVATE ees_workbench)

# Testing Infrastructure
enable_testing()

add_executable(test_core tests/test_core.cpp)
target_link_libraries(test_core PRIVATE ees_core)
add_test(NAME CoreMathTests COMMAND test_core)

add_executable(test_unit1 tests/test_unit1.cpp)
target_link_libraries(test_unit1 PRIVATE ees_core ees_unit1)
add_test(NAME Unit1Tests COMMAND test_unit1)

add_executable(test_unit2 tests/test_unit2.cpp)
target_link_libraries(test_unit2 PRIVATE ees_core ees_unit2)
add_test(NAME Unit2Tests COMMAND test_unit2)

add_executable(test_unit3 tests/test_unit3.cpp)
target_link_libraries(test_unit3 PRIVATE ees_core ees_unit3)
add_test(NAME Unit3Tests COMMAND test_unit3)

add_executable(test_unit4 tests/test_unit4.cpp)
target_link_libraries(test_unit4 PRIVATE ees_core ees_unit4)
add_test(NAME Unit4Tests COMMAND test_unit4)
