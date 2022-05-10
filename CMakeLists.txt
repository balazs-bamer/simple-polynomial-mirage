cmake_minimum_required(VERSION 3.16)
set(CMAKE_CXX_COMPILER "/usr/bin/g++")
project(RungeKuttaRayBending)

#enable_testing()
#find_package(GTest REQUIRED)

#add_compile_options(-ggdb -D_GLIBCXX_DEBUG -std=c++17)
add_compile_options(-O2 -std=c++17)

add_definitions(-DEIGEN_MATRIX_PLUGIN="Matrix_initializer_list.h" -DEIGEN_ARRAY_PLUGIN="Array_initializer_list.h")

INCLUDE_DIRECTORIES ( "/usr/include/eigen3"
                      "../repos/stl_reader"
                      "../repos/eigen-initializer_list/src" )
ADD_LIBRARY (RungeKuttaRayBendingLib SHARED RungeKuttaRayBending.cpp mathUtil.cpp)
target_link_libraries(RungeKuttaRayBendingLib quadmath png gsl pthread)

#add_executable(googleTest googleTest.cpp)
#target_link_libraries(googleTest GTest::GTest GTest::Main RungeKuttaRayBendingLib quadmath)
#add_test(google-test googleTest)

add_executable(main main.cpp simpleRaytracer.cpp)
target_link_libraries(main RungeKuttaRayBendingLib quadmath png gsl)
