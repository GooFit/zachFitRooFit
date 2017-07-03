cmake_minimum_required(VERSION 3.4)

project(ZachFitRooFit CXX)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" "${CMAKE_MODULE_PATH}")

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_package(ROOT 6 REQUIRED)

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${ROOTSYS}/cmake/modules")

if(NOT "${CMAKE_CURRENT_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_BINARY_DIR}")
    execute_process(COMMAND
        ${CMAKE_COMMAND} -E create_symlink
        "${CMAKE_CURRENT_SOURCE_DIR}/DstarWidthAnalysis_Data"
        "${CMAKE_CURRENT_BINARY_DIR}/DstarWidthAnalysis_Data"
        )
endif()


include_directories(.)

root_add_dictionary(
    RooRBWGaussConv_Dict
    RooRBWGaussConv.h)

root_add_dictionary(
    RooArgusGenBG_Dict
    RooArgusGenBG.h)


add_library(RooArgusGenBG SHARED RooArgusGenBG.cxx)
target_link_libraries(RooArgusGenBG PUBLIC ROOT::ROOT ROOT::RooFit ROOT::RooFitCore RooArgusGenBG_Dict)

add_library(RooRBWGaussConv SHARED RooRBWGaussConv.cxx)
target_link_libraries(RooRBWGaussConv PUBLIC ROOT::ROOT ROOT::RooFit ROOT::RooFitCore RooRBWGaussConv_Dict)

add_executable(zachFit_roofit zachFit_roofit.C)
target_link_libraries(zachFit_roofit PUBLIC RooArgusGenBG RooRBWGaussConv)
