/*
 * Catch C++ standard exceptions and raise corresponding
 * Python exceptions. Use simplified GIL API to ensure 
 * we have the GIL.
 */


#ifndef CKDTREE_CPP_EXC_H
#define CKDTREE_CPP_EXC_H

extern "C" void 
translate_cpp_exception(); 

extern "C" void 
translate_cpp_exception_with_gil();

#endif

