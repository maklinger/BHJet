cmake_minimum_required(VERSION 3.12)
project(PyBHJet)

# Specify CMAKE_PREFIX_PATH for pybind11
set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};/Users/ruby/anaconda3/envs/threeML/lib/python3.9/site-packages/pybind11/share/cmake/pybind11")

# Find pybind11
find_package(pybind11 REQUIRED)
find_package(GSL REQUIRED)
target_link_libraries(pybhjet PRIVATE GSL::GSL)
#add gsl include directories link 

# Add the library
add_library(pybhjet MODULE python_interface.cpp bhjet_class.cpp)

# Remove the 'lib' prefix and set the output name to match the module name
set_target_properties(pybhjet PROPERTIES PREFIX "")
set_target_properties(pybhjet PROPERTIES OUTPUT_NAME "pybhjet")

# Link pybind11 and add its include directories
target_link_libraries(pybhjet PRIVATE pybind11::module)
target_include_directories(pybhjet PRIVATE /Users/ruby/anaconda3/envs/threeML/lib/python3.9/site-packages/pybind11/include)