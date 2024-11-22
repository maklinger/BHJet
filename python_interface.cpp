#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/buffer_info.h>
#include <pybind11/functional.h>
#include "bhjet_class.h"

namespace py = pybind11; 

#define SET_ARGS_VEC(classtype, function, type) \
	[](classtype &a, py::array_t<type> array) { \
		py::buffer_info buf = array.request(); \
	    if (buf.ndim != 1) { \
	        throw std::runtime_error("Number of dimensions must be one and/ or array size is wrong"); \
		} \
		a.function(std::vector<type>(array.data(), array.data() + array.size()));}



PYBIND11_MODULE(pybhjet, pybhjetmodule){

    py::class_<BhJetClass> pbh(pybhjetmodule, "PyBHJet"); 
    pbh.def(py::init<double>()); 

    pbh.def("run", &BhJetClass::run, "Running BHJet C++");
    pbh.def("get_photon_array", &BhJetClass::get_photon_array);

    pbh.def("set_photon_array", SET_ARGS_VEC(BhJetClass, set_photon_array, double));

    
}