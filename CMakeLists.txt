cmake_minimum_required(VERSION 3.4)

project(ZachFitRooFit CXX)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" "${CMAKE_MODULE_PATH}")

find_package(ROOT 6 REQUIRED)

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${ROOTSYS}/cmake/modules")

include_directories(".")

function(root_add_dictionary OUTNAME)
    add_custom_command(
        OUTPUT   ${CMAKE_CURRENT_BINARY_DIR}/${OUTNAME}.cxx
        COMMAND ${ROOTCINT_EXECUTABLE} -v4 -f ${CMAKE_CURRENT_BINARY_DIR}/${OUTNAME}.cxx ${ARGN}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        DEPENDS ${ARN}
        )
    add_library(${OUTNAME} STATIC ${CMAKE_CURRENT_BINARY_DIR}/${OUTNAME}.cxx)
    target_link_libraries(${OUTNAME} PUBLIC ROOT::ROOT)
endfunction()

root_add_dictionary(
    RooRBWGaussConv_Dict
    RooRBWGaussConv.h)
#    Linkdef.h)

root_add_dictionary(
    RooArgusGenBG_Dict
    RooArgusGenBG.h)
#    Linkdef.h)


add_library(RooArgusGenBG SHARED RooArgusGenBG.cxx)
target_link_libraries(RooArgusGenBG PUBLIC ROOT::ROOT ROOT::RooFit ROOT::RooFitCore RooArgusGenBG_Dict)

add_library(RooRBWGaussConv SHARED RooRBWGaussConv.cxx)
target_link_libraries(RooRBWGaussConv PUBLIC ROOT::ROOT ROOT::RooFit ROOT::RooFitCore RooRBWGaussConv_Dict)

add_executable(zachFit_roofit zachFit_roofit.C)
target_link_libraries(zachFit_roofit PUBLIC RooRBWGaussConv RooRBWGaussConv)
