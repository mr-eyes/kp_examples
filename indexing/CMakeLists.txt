cmake_minimum_required (VERSION 3.4)
project (kTests C CXX)
set(default_build_type "Release")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++17 -fPIC -lgomp -lrt -fopenmp -W -Wall -pedantic -O3 -Ofast -march=native")


# ----------------------------------------------------------------------------
# kProcessor Setup
# ----------------------------------------------------------------------------

add_subdirectory(lib/kProcessor)
set(kProcessor_INCLUDE_PATH "${PROJECT_SOURCE_DIR}/lib/kProcessor/include")

# Add include directories.
include_directories("${PROJECT_SOURCE_DIR}/lib/kProcessor/ThirdParty/MQF/include")
include_directories("${PROJECT_SOURCE_DIR}/lib/kProcessor/ThirdParty/Blight")
include_directories (${kProcessor_INCLUDE_PATH})


# ----------------------------------------------------------------------------
# END kProcessor Setup
# ----------------------------------------------------------------------------


add_executable (test_indexing main.cpp)
target_link_libraries (test_indexing kProcessor z)
target_include_directories(test_indexing INTERFACE ${kProcessor_INCLUDE_PATH})

#add_executable (testKseq test_kseq.cpp)
#target_link_libraries (testKseq z)



# genesIndexing executable (Have an issue after some time in the execution)
#add_executable (genesIndexing genes_indexing.cpp)
#target_link_libraries (genesIndexing kProcessor pthread z)
#target_include_directories(genesIndexing INTERFACE ${kProcessor_INCLUDE_PATH})