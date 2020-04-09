import ctypes
import os

highslib = ctypes.cdll.LoadLibrary("libhighs.so") # highs lib folder must be in "LD_LIBRARY_PATH" environment variable

highslib.Highs_call.argtypes = (ctypes.c_int, ctypes.c_int, ctypes.c_int, 
ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), 
ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), 
ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double),
ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int))
highslib.Highs_call.restype = ctypes.c_int

def Highs_call(colcost, collower, colupper, rowlower, rowupper, astart, aindex, avalue):
   global highslib
   n_col = len(colcost)
   n_row = len(rowlower)
   n_nz = len(aindex)

   dbl_array_type_col = ctypes.c_double * n_col
   dbl_array_type_row = ctypes.c_double * n_row
   int_array_type_astart = ctypes.c_int * (n_col + 1)
   int_array_type_aindex = ctypes.c_int * n_nz
   dbl_array_type_avalue = ctypes.c_double * n_nz

   int_array_type_col = ctypes.c_int * n_col
   int_array_type_row = ctypes.c_int * n_row

   col_value = [0] * n_col
   col_dual = [0] * n_col

   row_value = [0] * n_col
   row_dual = [0] * n_col

   col_basis = [0] * n_col
   row_basis = [0] * n_row

   return_val = 0

   col_value = dbl_array_type_col(*col_value)
   col_dual = dbl_array_type_col(*col_dual)
   row_value = dbl_array_type_row(*row_value)
   row_dual = dbl_array_type_row(*row_dual)
   col_basis = int_array_type_col(*col_basis)
   row_basis = int_array_type_row(*row_basis)

   retcode = highslib.Highs_call(
      ctypes.c_int(n_col), ctypes.c_int(n_row), ctypes.c_int(n_nz), 
      dbl_array_type_col(*colcost), dbl_array_type_col(*collower), dbl_array_type_col(*colupper), 
      dbl_array_type_row(*rowlower), dbl_array_type_row(*rowupper), 
      int_array_type_astart(*astart), int_array_type_aindex(*aindex), dbl_array_type_avalue(*avalue),
      col_value, col_dual, 
      row_value, row_dual, 
      col_basis, row_basis, ctypes.byref(ctypes.c_int(return_val)))
   return retcode, list(col_value), list(col_dual), list(row_value), list(row_dual), list(col_basis), list(row_basis)

