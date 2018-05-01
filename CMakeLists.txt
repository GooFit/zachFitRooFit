cmake_minimum_required(VERSION 3.4)

project(ZachFitRooFit LANGUAGES CXX)

find_package(ROOT 6 CONFIG REQUIRED COMPONENTS RooFit)

include("${ROOT_DIR}/modules/RootNewMacros.cmake")

if(NOT "${CMAKE_CURRENT_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_BINARY_DIR}")
    execute_process(COMMAND
        ${CMAKE_COMMAND} -E create_symlink
        "${CMAKE_CURRENT_SOURCE_DIR}/DstarWidthAnalysis_Data"
        "${CMAKE_CURRENT_BINARY_DIR}/DstarWidthAnalysis_Data")
endif()

# Needed to be set on the directory in this way
# for root_generate_dictionary
include_directories("${CMAKE_CURRENT_SOURCE_DIR}")

root_generate_dictionary(
    G__RooRBWGaussConv
    RooRBWGaussConv.h
  LINKDEF
    RooRBWGaussConvLinkDef.h)

root_generate_dictionary(
    G__RooArgusGenBG
    RooArgusGenBG.h
  LINKDEF
    RooArgusGenBGLinkDef.h)

add_library(RooArgusGenBG SHARED RooArgusGenBG.cxx G__RooArgusGenBG.cxx)
target_link_libraries(RooArgusGenBG PUBLIC ROOT::Libraries ROOT::Flags)

add_library(RooRBWGaussConv SHARED RooRBWGaussConv.cxx G__RooRBWGaussConv.cxx)
target_link_libraries(RooRBWGaussConv PUBLIC ROOT::Libraries ROOT::Flags)

add_executable(zachFit_roofit zachFit_roofit.C)
target_link_libraries(zachFit_roofit PUBLIC RooArgusGenBG RooRBWGaussConv)
