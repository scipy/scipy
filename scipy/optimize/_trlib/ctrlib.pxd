# MIT License
#
# Copyright (c) 2016--2017 Felix Lenders
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

cimport libc.stdio

cdef extern from "trlib.h":
    cdef long _TRLIB_CLR_CONV_BOUND    "TRLIB_CLR_CONV_BOUND"    
    cdef long _TRLIB_CLR_CONV_INTERIOR "TRLIB_CLR_CONV_INTERIOR" 
    cdef long _TRLIB_CLR_APPROX_HARD   "TRLIB_CLR_APPROX_HARD"   
    cdef long _TRLIB_CLR_NEWTON_BREAK  "TRLIB_CLR_NEWTON_BREAK"  
    cdef long _TRLIB_CLR_HARD_INIT_LAM "TRLIB_CLR_HARD_INIT_LAM" 
    cdef long _TRLIB_CLR_CONTINUE      "TRLIB_CLR_CONTINUE"      
    cdef long _TRLIB_CLR_ITMAX         "TRLIB_CLR_ITMAX"         
    cdef long _TRLIB_CLR_FAIL_FACTOR   "TRLIB_CLR_FAIL_FACTOR"   
    cdef long _TRLIB_CLR_FAIL_LINSOLVE "TRLIB_CLR_FAIL_LINSOLVE" 
    cdef long _TRLIB_CLR_FAIL_TTR      "TRLIB_CLR_FAIL_TTR"      
    cdef long _TRLIB_CLR_PCINDEF       "TRLIB_CLR_PCINDEF"       
    cdef long _TRLIB_CLR_UNEXPECT_INT  "TRLIB_CLR_UNEXPECT_INT"  
    cdef long _TRLIB_CLR_FAIL_HARD     "TRLIB_CLR_FAIL_HARD"    
    cdef long _TRLIB_CLT_CG            "TRLIB_CLT_CG"            
    cdef long _TRLIB_CLT_L             "TRLIB_CLT_L"             
    cdef long _TRLIB_CLA_TRIVIAL       "TRLIB_CLA_TRIVIAL"       
    cdef long _TRLIB_CLA_INIT          "TRLIB_CLA_INIT"          
    cdef long _TRLIB_CLA_RETRANSF      "TRLIB_CLA_RETRANSF"      
    cdef long _TRLIB_CLA_UPDATE_STATIO "TRLIB_CLA_UPDATE_STATIO" 
    cdef long _TRLIB_CLA_UPDATE_GRAD   "TRLIB_CLA_UPDATE_GRAD"   
    cdef long _TRLIB_CLA_UPDATE_DIR    "TRLIB_CLA_UPDATE_DIR"    
    cdef long _TRLIB_CLA_NEW_KRYLOV    "TRLIB_CLA_NEW_KRYLOV"    
    cdef long _TRLIB_CLA_CONV_HARD     "TRLIB_CLA_CONV_HARD"     
    cdef long _TRLIB_CLA_OBJVAL        "TRLIB_CLA_OBJVAL"        
    cdef long _TRLIB_CLS_INIT          "TRLIB_CLS_INIT"          
    cdef long _TRLIB_CLS_HOTSTART      "TRLIB_CLS_HOTSTART"      
    cdef long _TRLIB_CLS_HOTSTART_G    "TRLIB_CLS_HOTSTART_G"    
    cdef long _TRLIB_CLS_HOTSTART_R    "TRLIB_CLS_HOTSTART_R"    
    cdef long _TRLIB_CLS_HOTSTART_T    "TRLIB_CLS_HOTSTART_T"    
    cdef long _TRLIB_CLS_HOTSTART_P    "TRLIB_CLS_HOTSTART_P"    
    cdef long _TRLIB_CLS_VEC_INIT      "TRLIB_CLS_VEC_INIT"      
    cdef long _TRLIB_CLS_CG_NEW_ITER   "TRLIB_CLS_CG_NEW_ITER"   
    cdef long _TRLIB_CLS_CG_UPDATE_S   "TRLIB_CLS_CG_UPDATE_S"   
    cdef long _TRLIB_CLS_CG_UPDATE_GV  "TRLIB_CLS_CG_UPDATE_GV"  
    cdef long _TRLIB_CLS_CG_UPDATE_P   "TRLIB_CLS_CG_UPDATE_P"   
    cdef long _TRLIB_CLS_LANCZOS_SWT   "TRLIB_CLS_LANCZOS_SWT"   
    cdef long _TRLIB_CLS_L_UPDATE_P    "TRLIB_CLS_L_UPDATE_P"    
    cdef long _TRLIB_CLS_L_CMP_CONV    "TRLIB_CLS_L_CMP_CONV"    
    cdef long _TRLIB_CLS_L_CMP_CONV_RT "TRLIB_CLS_L_CMP_CONV_RT" 
    cdef long _TRLIB_CLS_L_CHK_CONV    "TRLIB_CLS_L_CHK_CONV"    
    cdef long _TRLIB_CLS_L_NEW_ITER    "TRLIB_CLS_L_NEW_ITER"    
    cdef long _TRLIB_CLS_CG_IF_IRBLK_P "TRLIB_CLS_CG_IF_IRBLK_P" 
    cdef long _TRLIB_CLS_CG_IF_IRBLK_C "TRLIB_CLS_CG_IF_IRBLK_C" 
    cdef long _TRLIB_CLS_CG_IF_IRBLK_N "TRLIB_CLS_CG_IF_IRBLK_N" 
    cdef long _TRLIB_CLC_NO_EXP_INV    "TRLIB_CLC_NO_EXP_INV"    
    cdef long _TRLIB_CLC_EXP_INV_LOC   "TRLIB_CLC_EXP_INV_LOC"   
    cdef long _TRLIB_CLC_EXP_INV_GLO   "TRLIB_CLC_EXP_INV_GLO"   
    cdef long _TRLIB_CLT_CG_INT        "TRLIB_CLT_CG_INT"        
    cdef long _TRLIB_CLT_CG_BOUND      "TRLIB_CLT_CG_BOUND"      
    cdef long _TRLIB_CLT_LANCZOS       "TRLIB_CLT_LANCZOS"       
    cdef long _TRLIB_CLT_HOTSTART      "TRLIB_CLT_HOTSTART"      

    long trlib_krylov_prepare_memory(long itmax, double *fwork)
    long trlib_krylov_memory_size(long itmax, long *iwork_size, long *fwork_size, long *h_pointer)
    long trlib_krylov_min(
        long init, double radius, long equality, long itmax, long itmax_lanczos,
        double tol_rel_i, double tol_abs_i,
        double tol_rel_b, double tol_abs_b, double zero, double obj_lo,
        long ctl_invariant, long convexify, long earlyterm,
        double g_dot_g, double v_dot_g, double p_dot_Hp,
        long *iwork, double *fwork, int refine,
        long verbose, long unicode, char *prefix, libc.stdio.FILE *fout, long *timing,
        long *action, long *iter, long *ityp,
        double *flt1, double *flt2, double *flt3)
    long trlib_krylov_timing_size()
