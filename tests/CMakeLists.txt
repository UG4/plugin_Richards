# Richards unit tests
cmake_minimum_required(VERSION 3.16)
project(RichardsTests LANGUAGES CXX)

enable_testing()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(GTest REQUIRED)

add_executable(richards_test
    ../van_genuchten.cpp
    richards_unit_test.cpp
)

target_include_directories(richards_test PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../../..
    ${CMAKE_CURRENT_SOURCE_DIR}/../../../externals/JSONForUG4/json-cxx/single_include
    ${CMAKE_CURRENT_SOURCE_DIR}/../../../externals/AutodiffForUG4/autodiff
    ${CMAKE_CURRENT_SOURCE_DIR}/../../../ugcore/ugbase
)

target_link_directories(richards_test
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../../../lib)

target_link_libraries(richards_test 
    PRIVATE GTest::gtest_main ug4)

add_test(NAME RichardsUnitTests COMMAND richards_test)
