cmake_minimum_required(VERSION 3.15...3.30)

project(${SKBUILD_PROJECT_NAME} LANGUAGES CXX)

find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)

Python_add_library(_uarray MODULE
    small_dynamic_array.h
    vectorcall.h
    vectorcall.cxx
    _uarray_dispatch.cxx
WITH_SOABI)
set_property(TARGET _uarray PROPERTY CXX_STANDARD 17)
install(TARGETS _uarray LIBRARY DESTINATION uarray)
