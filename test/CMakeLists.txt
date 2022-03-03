cmake_minimum_required(VERSION 3.18)
project(GQTest VERSION 1.0.0)

# name your project, set executable and link SingGQ to your executable
add_executable(GQTest main.cpp)
find_package(SingGQ REQUIRED)
target_link_libraries(GQTest PUBLIC SingGQ::SingGQ)
