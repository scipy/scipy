#ifndef COMPLEX_OPS_H
#define COMPLEX_OPS_H

/*
 *  Functions to handle arithmetic operations on NumPy complex values
 */




/*
 * Addition
 */
inline npy_cfloat operator+(const npy_cfloat& A, const npy_cfloat& B){
  npy_cfloat result;
  result.real = A.real + B.real;
  result.imag = A.imag + B.imag;
  return result;
}
inline npy_cdouble operator+(const npy_cdouble& A, const npy_cdouble& B){
  npy_cdouble result;
  result.real = A.real + B.real;
  result.imag = A.imag + B.imag;
  return result;
}
inline npy_clongdouble operator+(const npy_clongdouble& A, const npy_clongdouble& B){
  npy_clongdouble result;
  result.real = A.real + B.real;
  result.imag = A.imag + B.imag;
  return result;
}

inline npy_cfloat& operator+=(npy_cfloat& A, const npy_cfloat& B){
  A.real += B.real;
  A.imag += B.imag;
  return A;
}
inline npy_cdouble& operator+=(npy_cdouble& A, const npy_cdouble& B){
  A.real += B.real;
  A.imag += B.imag;
  return A;
}
inline npy_clongdouble& operator+=(npy_clongdouble& A, const npy_clongdouble& B){
  A.real += B.real;
  A.imag += B.imag;
  return A;
}



/*
 * Subtraction
 */
inline npy_cfloat operator-(const npy_cfloat& A, const npy_cfloat& B){
  npy_cfloat result;
  result.real = A.real - B.real;
  result.imag = A.imag - B.imag;
  return result;
}
inline npy_cdouble operator-(const npy_cdouble& A, const npy_cdouble& B){
  npy_cdouble result;
  result.real = A.real - B.real;
  result.imag = A.imag - B.imag;
  return result;
}
inline npy_clongdouble operator-(const npy_clongdouble& A, const npy_clongdouble& B){
  npy_clongdouble result;
  result.real = A.real - B.real;
  result.imag = A.imag - B.imag;
  return result;
}

inline npy_cfloat& operator-=(npy_cfloat& A, const npy_cfloat& B){
  A.real -= B.real;
  A.imag -= B.imag;
  return A;
}
inline npy_cdouble& operator-=(npy_cdouble& A, const npy_cdouble& B){
  A.real -= B.real;
  A.imag -= B.imag;
  return A;
}
inline npy_clongdouble& operator-=(npy_clongdouble& A, const npy_clongdouble& B){
  A.real -= B.real;
  A.imag -= B.imag;
  return A;
}


/*
 * Multiplication
 */
inline npy_cfloat operator*(const npy_cfloat& A, const npy_cfloat& B){
  npy_cfloat result;
  result.real = A.real * B.real - A.imag * B.imag;
  result.imag = A.real * B.imag + A.imag * B.real;
  return result;
}
inline npy_cdouble operator*(const npy_cdouble& A, const npy_cdouble& B){
  npy_cdouble result;
  result.real = A.real * B.real - A.imag * B.imag;
  result.imag = A.real * B.imag + A.imag * B.real;
  return result;
}
inline npy_clongdouble operator*(const npy_clongdouble& A, const npy_clongdouble& B){
  npy_clongdouble result;
  result.real = A.real * B.real - A.imag * B.imag;
  result.imag = A.real * B.imag + A.imag * B.real;
  return result;
}

inline npy_cfloat& operator*=(npy_cfloat& A, const npy_cfloat& B){
  npy_float temp = A.real * B.real - A.imag * B.imag;
  A.imag = A.real * B.imag + A.imag * B.real;
  A.real = temp;
  return A;
}
inline npy_cdouble& operator*=(npy_cdouble& A, const npy_cdouble& B){
  npy_double temp = A.real * B.real - A.imag * B.imag;
  A.imag = A.real * B.imag + A.imag * B.real;
  A.real = temp;
  return A;
}
inline npy_clongdouble& operator*=(npy_clongdouble& A, const npy_clongdouble& B){
  npy_longdouble temp = A.real * B.real - A.imag * B.imag;
  A.imag = A.real * B.imag + A.imag * B.real;
  A.real = temp;
  return A;
}


/*
 * Division
 */
inline npy_cfloat operator/(const npy_cfloat& A, const npy_cfloat& B){
  npy_cfloat result;
  npy_float denom = 1.0 / (B.real * B.real + B.imag * B.imag);
  result.real = (A.real * B.real + A.imag * B.imag) * denom;
  result.imag = (A.real * B.imag - A.imag * B.real) * denom;
  return result;
}
inline npy_cdouble operator/(const npy_cdouble& A, const npy_cdouble& B){
  npy_cdouble result;
  npy_double denom = 1.0 / (B.real * B.real + B.imag * B.imag);
  result.real = (A.real * B.real + A.imag * B.imag) * denom;
  result.imag = (A.real * B.imag - A.imag * B.real) * denom;
  return result;
}
inline npy_clongdouble operator/(const npy_clongdouble& A, const npy_clongdouble& B){
  npy_clongdouble result;
  npy_longdouble denom = 1.0 / (B.real * B.real + B.imag * B.imag);
  result.real = (A.real * B.real + A.imag * B.imag) * denom;
  result.imag = (A.real * B.imag - A.imag * B.real) * denom;
  return result;
}

inline npy_cfloat& operator/=(npy_cfloat& A, const npy_cfloat& B){
  A = A*B;
  return A;
}
inline npy_cdouble& operator/=(npy_cdouble& A, const npy_cdouble& B){
  A = A*B;
  return A;
}
inline npy_clongdouble& operator/=(npy_clongdouble& A, const npy_clongdouble& B){
  A = A*B;
  return A;
}

/*
 * Equality (complex==complex)
 */
inline bool operator==(const npy_cfloat& A, const npy_cfloat& B){
  return A.real == B.real && A.imag == B.imag;
}
inline bool operator==(const npy_cdouble& A, const npy_cdouble& B){
  return A.real == B.real && A.imag == B.imag;
}
inline bool operator==(const npy_clongdouble& A, const npy_clongdouble& B){
  return A.real == B.real && A.imag == B.imag;
}

inline bool operator!=(const npy_cfloat& A, const npy_cfloat& B){
  return A.real != B.real || A.imag != B.imag;
}
inline bool operator!=(const npy_cdouble& A, const npy_cdouble& B){
  return A.real != B.real || A.imag != B.imag;
}
inline bool operator!=(const npy_clongdouble& A, const npy_clongdouble& B){
  return A.real != B.real || A.imag != B.imag;
}

/*
 * Equality (complex==scalar)
 */
inline bool operator==(const npy_cfloat& A, const npy_float& B){
  return A.real == B && A.imag == 0;
}
inline bool operator==(const npy_cdouble& A, const npy_double& B){
  return A.real == B && A.imag == 0;
}
inline bool operator==(const npy_clongdouble& A, const npy_longdouble& B){
  return A.real == B && A.imag == 0;
}

inline bool operator!=(const npy_cfloat& A, const npy_float& B){
  return A.real != B || A.imag != 0;
}
inline bool operator!=(const npy_cdouble& A, const npy_double& B){
  return A.real != B || A.imag != 0;
}
inline bool operator!=(const npy_clongdouble& A, const npy_longdouble& B){
  return A.real != B || A.imag != 0;
}

/*
 * Equality (scalar==complex)
 */
inline bool operator==(const npy_float& A, const npy_cfloat& B){
  return A == B.real && 0 == B.imag;
}
inline bool operator==(const npy_double& A, const npy_cdouble& B){
  return A == B.real && 0 == B.imag;
}
inline bool operator==(const npy_longdouble& A, const npy_clongdouble& B){
  return A == B.real && 0 == B.imag;
}

inline bool operator!=(const npy_float& A, const npy_cfloat& B){
  return A != B.real || 0 != B.imag;
}
inline bool operator!=(const npy_double& A, const npy_cdouble& B){
  return A != B.real || 0 != B.imag;
}
inline bool operator!=(const npy_longdouble& A, const npy_clongdouble& B){
  return A != B.real || 0 != B.imag;
}

#endif
