cmake_minimum_required(VERSION 3.16)
project(csr CXX)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

#库本体
add_library(csr src/matrix.cpp
                src/csr_matrix.cpp 
                src/poisson.cpp 
                src/linalg.cpp 
                src/grid.cpp
                src/flux.cpp
                src/advect1d.cpp
                src/burgers.cpp)
target_include_directories(csr PUBLIC include)

#---GoogleTest (FetchContent 方式) ----
include(FetchContent)
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.14.0   #锁定版本，保证可复现
) 
#对Windows有用，Linux 上无害，官方推荐写法
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

#可执行文件：薄薄的一层入口
add_executable(test_matrix tests/test_matrix.cpp)
target_link_libraries(test_matrix PRIVATE csr GTest::gtest_main)
add_test(NAME matrix_tests COMMAND test_matrix)

add_executable(test_csr_matrix tests/test_csr_matrix.cpp)
target_link_libraries(test_csr_matrix PRIVATE csr GTest::gtest_main)
add_test(NAME csr_matrix_tests COMMAND test_csr_matrix)

add_executable(test_poisson tests/test_poisson.cpp)
target_link_libraries(test_poisson PRIVATE csr GTest::gtest_main)
add_test(NAME poisson_tests COMMAND test_poisson)

add_executable(test_convergence tests/test_convergence.cpp)
target_link_libraries(test_convergence PRIVATE csr GTest::gtest_main)
add_test(NAME convergence_tests COMMAND test_convergence)

add_executable(test_grid tests/test_grid.cpp)
target_link_libraries(test_grid PRIVATE csr GTest::gtest_main)
add_test(NAME grid_tests COMMAND test_grid)

add_executable(test_flux tests/test_flux.cpp)
target_link_libraries(test_flux PRIVATE csr GTest::gtest_main)
add_test(NAME flux_tests COMMAND test_flux)

add_executable(test_advection tests/test_advection.cpp)
target_link_libraries(test_advection PRIVATE csr GTest::gtest_main)
add_test(NAME advection_tests COMMAND test_advection)

add_executable(test_burgers tests/test_burgers.cpp)
target_link_libraries(test_burgers PRIVATE csr GTest::gtest_main)
add_test(NAME burgers_tests COMMAND test_burgers)

add_executable(advect_demo drivers/advect_demo.cpp)
target_link_libraries(advect_demo PRIVATE csr)

include(GoogleTest)
gtest_discover_tests(test_matrix 
                     test_csr_matrix
                     test_poisson 
                     test_convergence
                     test_grid
                     test_flux
                     test_advection)