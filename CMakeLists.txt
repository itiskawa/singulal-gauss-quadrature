cmake_minimum_required(VERSION 3.2)

if( TARGET GQJacobi)
    return()
endif(TARGET GQJAcobi)

if( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR} )
  set(STANDALONE ON)
endif( )


project(GQJacobi)
message(STATUS "GQJacobi")

find_package(Eigen3 REQUIRED)

string( REGEX REPLACE "^lib" "" LIB_NAME ${PROJECT_NAME} )
add_library( ${LIB_NAME} INTERFACE )
add_library( ${PROJECT_NAME}::${LIB_NAME} ALIAS ${LIB_NAME} )
target_include_directories( ${LIB_NAME} INTERFACE
  $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/src>
  $<BUILD_INTERFACE:${${PROJECT_NAME}_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:include/${PROJECT_NAME}>
)

target_sources(
    ${LIB_NAME}
    INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/GQJacobi/gauss-jacobi.hpp>
)

target_link_libraries(${LIB_NAME} INTERFACE Eigen3::Eigen)
target_compile_features(${LIB_NAME} INTERFACE cxx_std_11 )
