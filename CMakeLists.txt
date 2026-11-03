cmake_minimum_required(VERSION 3.16)
project(csr CXX)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

#库本体
add_library(csr src/matrix.cpp src/csr_matrix.cpp)
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
add_executable(test_csr_matrix tests/test_csr_matrix.cpp)
target_link_libraries(test_matrix PRIVATE csr GTest::gtest_main)
target_link_libraries(test_csr_matrix PRIVATE csr GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(test_matrix test_csr_matrix)