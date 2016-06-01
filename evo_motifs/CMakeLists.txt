cmake_minimum_required (VERSION 2.6)
project (EVO_DNA)

SET (PLATFORM 64)

# The version number.
set (EVO_DNA_VERSION_MAJOR 1)
set (EVO_DNA_VERSION_MINOR 0)

SET (CMAKE_CXX_FLAGS                "-Wall -std=c++11 -O3 -funroll-loops -msse3")
SET (CMAKE_CXX_FLAGS_DEBUG          "-g")
SET (CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
SET (CMAKE_CXX_FLAGS_RELEASE        "-O4 -funrol-loops -DNDEBUG")

set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /opt/local/lib)

include_directories(${PROJECT_SOURCE_DIR})
#message("project source dir = ${PROJECT_SOURCE_DIR}")

#This will look for the undvc_common library, which should be in the
#same directory as this project. It should be compiled in the ../undvc_common/build
#directory
FIND_PATH(UNDVC_COMMON_INCLUDE_DIR file_io.hxx
    ../../undvc_common         #add addiitonal paths here if undvc_common is in a different directory
    )
MESSAGE(STATUS "UNDVC_COMMON include directory: ${UNDVC_COMMON_INCLUDE_DIR}")

FIND_LIBRARY(UNDVC_COMMON_LIBRARY
    NAMES undvc_common
    PATHS ../../undvc_common/build         #add additional paths here if undvc_common is in a different directory
                                        #or if it has been build in a different directory
    PATH_SUFFIXES lib
    )
MESSAGE(STATUS "UNDVC_COMMON library: ${UNDVC_COMMON_LIBRARY}")

include_directories(${UNDVC_COMMON_INCLUDE_DIR})

#This will look for the tao library, which should be in the
#same directory as this project. It should be compiled in the ../tao/build
#directory
FIND_PATH(TAO_INCLUDE_DIR ./examples/benchmarks.hxx
    ../../tao #add addiitonal paths here if undvc_common is in a different directory
    )
MESSAGE(STATUS "TAO include directory: ${TAO_INCLUDE_DIR}")

FIND_LIBRARY(TAO_MPI_LIBRARY NAMES mpi_algorithms PATHS ../../tao/build/mpi PATH_SUFFIXES lib)
FIND_LIBRARY(TAO_SYNCH_LIBRARY NAMES synchronous_algorithms PATHS ../../tao/build/synchronous_algorithms PATH_SUFFIXES lib)
FIND_LIBRARY(TAO_ASYNCH_LIBRARY NAMES asynchronous_algorithms PATHS ../../tao/build/asynchronous_algorithms PATH_SUFFIXES lib)
FIND_LIBRARY(TAO_NN_LIBRARY NAMES neural_networks PATHS ../../tao/build/neural_networks PATH_SUFFIXES lib)
FIND_LIBRARY(TAO_UTIL_LIBRARY NAMES tao_util PATHS ../../tao/build/util  PATH_SUFFIXES lib)

set(TAO_LIBRARY ${TAO_MPI_LIBRARY} ${TAO_SYNCH_LIBRARY} ${TAO_ASYNCH_LIBRARY} ${TAO_NN_LIBRARY} ${TAO_UTIL_LIBRARY})
MESSAGE(STATUS "TAO library: ${TAO_LIBRARY}")

include_directories(${TAO_INCLUDE_DIR})


#UNDVC_COMMON has the FindBOINC.cmake and FindMYSQL.cmake files
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${UNDVC_COMMON_INCLUDE_DIR}/cmake/Modules/")

set(BOOST_LIB_PREFIX "lib")
set(Boost_LIB_PREFIX "lib")
set(BOOST_COMPILER "-vc100")
set(Boost_COMPILER "-vc100")
set(Boost_ADDITIONAL_VERSIONS "1.47" "1.47.0")

set(USING_HODOR "OFF")
option(USING_HODOR "USING_HODOR")

message(STATUS "USING HODOR: ${USING_HODOR}")

if (USING_HODOR MATCHES "ON")
#need these for hodor
    message(STATUS "USING HODOR! : ${USING_HODOR}")
    set (MPI_COMPILER "/opt/mvapich2-x/gnu/bin/mpic++")
    set (MPI_INCLUDE_PATH "/opt/mvapich2-x/gnu/include")
endif ()


find_package(Boost REQUIRED system filesystem)
find_package(MPI)

include_directories(${Boost_INCLUDE_DIR})
include_directories(${MPI_INCLUDE_PATH})

add_executable(motif_fitness motif_fitness)
target_link_libraries(motif_fitness ${MPI_LIBRARY} ${MPI_EXTRA_LIBRARY} ${UNDVC_COMMON_LIBRARY} ${TAO_LIBRARY} ${Boost_LIBRARIES})

add_executable(brute_force brute_force)
target_link_libraries(brute_force ${UNDVC_COMMON_LIBRARY} ${Boost_LIBRARIES})

