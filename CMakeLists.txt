cmake_minimum_required(VERSION 3.16)
project(csr CXX)
set (CMAKE_CXX_STANDARD 17)

#库：可服用的核心代码
add_library(csr src/matrix.cpp)
target_include_directories(csr PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

#可执行文件：薄薄的一层入口
add_executable(test_matrix tests/test_matrix.cpp)
target_link_libraries(test_matrix PRIVATE csr)

enable_testing()
add_test(NAME matrix COMMAND test_matrix)