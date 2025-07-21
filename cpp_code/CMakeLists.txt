cmake_minimum_required(VERSION 3.12)
project(PyBHJet)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Specify CMAKE_PREFIX_PATH for pybind11 
find_package(pybind11 REQUIRED) #should find the CMAKE_PREFIX_PATH for pybind11 in environment, can set it directly if it is not found. 
# use this in terminal to find path:  python -m pybind11 --cmakedir

# set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};/path/to/pybind11/share/cmake/pybind11")
find_package(GSL REQUIRED)

# environment variable for Python paths
if(NOT DEFINED ENV{PYTHONPATH})
    message(WARNING "PYTHONPATH environment variable is not set. Ensure pybind11 is accessible.")
endif()



# Fetch Kariba from your fork
include(FetchContent)
FetchContent_Declare(
  kariba
  GIT_REPOSITORY https://github.com/maklinger/kariba.git
  GIT_TAG main
)
FetchContent_MakeAvailable(kariba)


# Collect Kariba sources from your fork
file(GLOB KARIBA_SOURCES "${kariba_SOURCE_DIR}/src/*.cpp")
include_directories(${kariba_SOURCE_DIR}/src)
include_directories(${kariba_SOURCE_DIR}/src/kariba)


file(GLOB BHJET_SOURCES "${kariba_SOURCE_DIR}/examples/model/*.cpp")
include_directories(${kariba_SOURCE_DIR}/examples/model)


# # Add Kariba source files
# set(KARIBA_SOURCES
#     ${CMAKE_CURRENT_SOURCE_DIR}/Kariba/BBody.cpp
#     ${CMAKE_CURRENT_SOURCE_DIR}/Kariba/Cyclosyn.cpp
#     ${CMAKE_CURRENT_SOURCE_DIR}/Kariba/Kappa.cpp
#     ${CMAKE_CURRENT_SOURCE_DIR}/Kariba/Mixed.cpp
#     ${CMAKE_CURRENT_SOURCE_DIR}/Kariba/Particles.cpp
#     ${CMAKE_CURRENT_SOURCE_DIR}/Kariba/ShSDisk.cpp
#     ${CMAKE_CURRENT_SOURCE_DIR}/Kariba/Compton.cpp
#     ${CMAKE_CURRENT_SOURCE_DIR}/Kariba/Bknpower.cpp
#     ${CMAKE_CURRENT_SOURCE_DIR}/Kariba/Powerlaw.cpp
#     ${CMAKE_CURRENT_SOURCE_DIR}/Kariba/Radiation.cpp
#     ${CMAKE_CURRENT_SOURCE_DIR}/Kariba/Thermal.cpp
# )

# Add bhjet and pybind source files
set(PYBHJET_SOURCES
    python_interface.cpp
    # bhjet_class.cpp
    # jetmain.cpp
    # utils.cpp
    # jetpars.cpp
)

# merge all src files into one name
set(ALL_SOURCES ${PYBHJET_SOURCES} ${KARIBA_SOURCES} ${BHJET_SOURCES})

# Add pybind11 module
pybind11_add_module(pybhjet MODULE ${ALL_SOURCES})

# Link libraries
target_link_libraries(pybhjet PRIVATE GSL::gsl GSL::gslcblas m pybind11::module)

# Include directories
target_include_directories(pybhjet PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR} # for headers in cpp_code
    ${kariba_SOURCE_DIR}/src # for kariba headers
    ${kariba_SOURCE_DIR}/src/kariba # for kariba headers
    ${kariba_SOURCE_DIR}/examples/model # for kariba headers
    /opt/local/include
)

# Set library search paths
target_link_directories(pybhjet PRIVATE /opt/local/lib)

# set_target_properties(pybhjet PROPERTIES
#     LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/src/pybhjet"
# )

install(TARGETS pybhjet DESTINATION pybhjet)
