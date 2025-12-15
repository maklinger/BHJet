cmake_minimum_required(VERSION 3.15)
project(BHJet LANGUAGES CXX)

# --- Options ---
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# I added this only for unsuccesful profiling attempts
# set(CMAKE_BUILD_TYPE Debug)
# add_compile_options(-g)

# --- Dependencies ---
set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)
find_package(GSL REQUIRED)

# Fetch Kariba from github
include(FetchContent)
FetchContent_Declare(
  kariba
  GIT_REPOSITORY https://github.com/maklinger/kariba.git
  GIT_TAG bhjetnew
#   GIT_REPOSITORY https://github.com/antonpannekoek/kariba.git
#   GIT_TAG main
)
FetchContent_MakeAvailable(kariba)

# Collect Kariba sources from github
file(GLOB KARIBA_SOURCES "${kariba_SOURCE_DIR}/src/*.cpp")
include_directories(${kariba_SOURCE_DIR}/src)


add_library(kariba_lib STATIC ${KARIBA_SOURCES})
target_include_directories(kariba_lib PUBLIC
    ${kariba_SOURCE_DIR}/src/kariba
)


# --- Sources ---
set(SRC_FILES
    src/JetDynamics.cpp
    src/BLJet.cpp
    src/RadiationZone.cpp
    src/BHJet.cpp
    src/utils.cpp
)

# --- Python module ---
pybind11_add_module(bhjet
    libpython/bindings.cpp
    ${SRC_FILES}
)


# Include bhjet headers
target_include_directories(bhjet PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Link Kariba and other libraries
target_link_libraries(bhjet PRIVATE kariba_lib GSL::gsl GSL::gslcblas pybind11::module)
if(UNIX)
    target_link_libraries(bhjet PRIVATE m)
endif()


# --- Output ---
# installs library to /python/bhjet
set_target_properties(bhjet PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python/bhjet
)


install(TARGETS bhjet DESTINATION python)

