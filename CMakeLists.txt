cmake_minimum_required(VERSION 3.16)
set(CMAKE_CXX_COMPILER "/usr/bin/clang++")
project(ShepardRayBending)

enable_testing()
find_package(GTest REQUIRED)

#add_compile_options(-ggdb -D_GLIBCXX_DEBUG -std=c++17)
add_compile_options(-O2 -std=c++17)

add_definitions(-DEIGEN_MATRIX_PLUGIN="Matrix_initializer_list.h" -DEIGEN_ARRAY_PLUGIN="Array_initializer_list.h")

INCLUDE_DIRECTORIES ( "/usr/include/eigen3"
                      "../repos/stl_reader"
                      "../repos/eigen-initializer_list/src" )
ADD_LIBRARY (ShepardRayBendingLib SHARED ShepardRayBending.cpp mathUtil.cpp)

add_executable(googleTest googleTest.cpp)
target_link_libraries(googleTest GTest::GTest GTest::Main ShepardRayBendingLib quadmath)
add_test(google-test googleTest)

add_executable(main main.cpp simpleRaytracer.cpp)
target_link_libraries(main ShepardRayBendingLib quadmath png)
