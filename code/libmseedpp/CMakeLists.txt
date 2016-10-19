cmake_minimum_required(VERSION 2.8.3)
project(libmseed)

## Sources
set(libmseed_SRCS src/fileutils.c
                  src/genutils.c
                  src/gswap.c
                  src/lmplatform.c
                  src/lookup.c
                  src/msrutils.c
                  src/pack.c
                  src/packdata.c
                  src/traceutils.c
                  src/tracelist.c
                  src/parseutils.c
                  src/unpack.c
                  src/unpackdata.c
                  src/selection.c
                  src/logging.c
)

set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -std=c++11)
## Include headers
include_directories(inc)

## Add libmseed library
add_library(${PROJECT_NAME} ${libmseed_SRCS})

# test
add_executable(parse_one_file test/parse_one_file.c)
add_dependencies(parse_one_file ${PROJECT_NAME})
target_link_libraries(parse_one_file ${PROJECT_NAME})


