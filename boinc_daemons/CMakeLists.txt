cmake_minimum_required (VERSION 2.6)
project (DNA_AT_HOME)
# The version number.
set (DNA_AT_HOME_VERSION_MAJOR 1)
set (DNA_AT_HOME_VERSION_MINOR 0)

# These will make cmake use clang and clang++ instead of gcc and g++
#SET (CMAKE_C_COMPILER             "/usr/bin/clang")
SET (CMAKE_C_FLAGS                "-Wall -g -O2 -std=c99 -Wl,-Bstatic")
SET (CMAKE_CXX_FLAGS              "-Wall -g -O2 -std=c++0x")

unset(Boost_INCLUDE_DIR CACHE)
unset(Boost_LIBRARY_DIRS CACHE)

#This may need to be updated or Linux/Windows
set(CMAKE_LIBRARY_PATH /home/tdesell/boost_1_53_0/stage/lib ${CMAKE_LIBRARY_PATH} /Users/deselt/Dropbox/software/opencv/lib /opt/local/lib)

include_directories (${PROJECT_SOURCE_DIR})

#This will look for the undvc_common library, which should be in the
#same directory as TAO. It should be compiled in the ../undvc_common/build
#directory
FIND_PATH(UNDVC_COMMON_INCLUDE_DIR file_io.hxx
    ~/undvc_common         #add addiitonal paths here if undvc_common is in a different directory
    )
MESSAGE(STATUS "UNDVC_COMMON include directory: ${UNDVC_COMMON_INCLUDE_DIR}")

FIND_LIBRARY(UNDVC_COMMON_LIBRARY
    NAMES undvc_common
    PATHS ~/undvc_common/build         #add additional paths here if undvc_common is in a different directory
                                        #or if it has been build in a different directory
    PATH_SUFFIXES lib
    )
MESSAGE(STATUS "UNDVC_COMMON library: ${UNDVC_COMMON_LIBRARY}")

include_directories(${UNDVC_COMMON_INCLUDE_DIR})

#UNDVC_COMMON has the FindBOINC.cmake and FindMYSQL.cmake files
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${UNDVC_COMMON_INCLUDE_DIR}/cmake/Modules/")

find_package(Boost REQUIRED COMPONENTS filesystem system)
message(STATUS ${Boost_LIBRARIES})

find_package(MySQL REQUIRED)
find_package(BOINC REQUIRED)

find_path(BOINC_INCLUDE_DIR config.h PATHS /boinc/src/boinc/ /boinc/)

add_definitions( -D_BOINC_ )

find_package(OpenSSL REQUIRED)

include_directories(
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/../include
    ${PROJECT_SOURCE_DIR}/undvc_common
	${BOINC_INCLUDE_DIR}
    ${BOINC_INCLUDE_DIR}/api
    ${BOINC_INCLUDE_DIR}/db
    ${BOINC_INCLUDE_DIR}/lib
    ${BOINC_INCLUDE_DIR}/sched
	${MYSQL_INCLUDE_DIR}
    ${Boost_INCLUDE_DIR}
    ${BOINC_INCLUDE_DIR}/tools/
)


add_executable(dna_validator
        ${BOINC_INCLUDE_DIR}/sched/validator
        ${BOINC_INCLUDE_DIR}/sched/validate_util
        ${BOINC_INCLUDE_DIR}/sched/validate_util2
        dna_validation_policy)
target_link_libraries(dna_validator ${UNDVC_COMMON_LIBRARY} ${BOINC_SERVER_LIBRARIES} ${MYSQL_LIBRARIES})

add_executable(dna_assimilator
   ${BOINC_INCLUDE_DIR}/sched/assimilator
   ${BOINC_INCLUDE_DIR}/sched/validate_util
   dna_assimilation_policy)
target_link_libraries(dna_assimilator ${UNDVC_COMMON_LIBRARY} ${BOINC_SERVER_LIBRARIES} ${MYSQL_LIBRARIES} ${Boost_LIBRARIES})

add_executable(dna_work_generator dna_work_generator)
target_link_libraries(dna_work_generator ${UNDVC_COMMON_LIBRARY}
    ${BOINC_SERVER_LIBRARIES} ${MYSQL_LIBRARIES} ${OPENSSL_LIBRARIES} ${Boost_LIBRARIES})
