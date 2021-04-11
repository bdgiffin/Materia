#ifndef MATRIX_FUNCTION_H
#define MATRIX_FUNCTION_H

#include <math.h>

// evaluate the matrix-valued function of the incoming 2x2 symmetric matrix,
// and return via the incoming components
#pragma acc routine seq
template<typename Functor>
inline void matrix_function(float& a11, float& a22, float& a12, Functor f) {
  // compute the difference in the eigenvalues [5 flops, 1 sqrt]
  float gap = sqrtf((a11-a22)*(a11-a22)+4.0f*a12*a12);

  // compute the smallest eigenvalue [2 flops]
  float lam1 = 0.5f*(a11+a22-gap);

  // compute f() of the smallest eigenvalue [1 f()]
  float f1 = f(lam1);

  // compute the stabilized finite difference of f() [3 flops, 1 f(), 1 div, 1 fmax]
  float diff = (f(lam1+gap)-f1)/fmax(gap,lam1*FLT_EPSILON);

  // compute the coefficient for the identity term [2 flops]
  float f0 = f1-lam1*diff;

  // return the matrix function f(A) = f0*I + diff*A [5 flops]
  a11 = diff*a11+f0;
  a22 = diff*a22+f0;
  a12 *= diff;

  // operation totals:
  // 17 flops
  // 1 sqrtf()
  // 1 fmax()
  // 1 div
  // 2 f()
  
} // matrix_function()

// evaluate the tensor-valued jacobian of the incoming 2x2 symmetric matrix
// and return the matrix-valued function via the incoming components
#pragma acc routine seq
template<typename Functor>
inline void matrix_jacobian(float& a11, float& a22, float& a12,
			    float& j11, float& j22, float& j33,
			    float& j23, float& j13, float& j12,
			    Functor f, Functor dfdx) {
  // compute the difference in the eigenvalues [5 flops, 1 sqrt]
  float gap = sqrtf((a11-a22)*(a11-a22)+4.0f*a12*a12);

  // compute the smallest eigenvalue [2 flops]
  float lam1 = 0.5f*(a11+a22-gap);

  // compute the stabilized inverse of the eigenvalue gap [1 flop, 1 div, 1 fmax]
  float igap = 1.0/fmax(gap,lam1*FLT_EPSILON);

  // compute f() of the smallest eigenvalue [1 f()]
  float f1 = f(lam1);

  // compute the stabilized finite difference of f() [2 flops, 1 f()]
  float lam2 = lam1+gap
  float diff = f(lam2)-f1;

  // compute projection matrix P2 = (A - lam1*I)*igap [5 flops]
  a11 = (a11-lam1)*igap;
  a22 = (a22-lam1)*igap;
  a12 *= igap;

  // compute jacobian coefficients [7 flops, 2 dfdx()]
  float j0 = dfdx(lam1);
  float j1 = ((j0+dfdx(lam2))*gap-2.0f*diff)*igap;
  float j2 = (diff-j0*gap)*igap;

  // compute the jacobian entries [20 flops]
  j11 = j0+a11*(j1*a11+2.0f*j2);
  j22 = j0+a22*(j1*a22+2.0f*j2);
  j12 = j1*a12*a12;
  j33 = j0+j12+j2*(a11+a22);
  j13 = a12*(j1*a11+j2);
  j23 = a12*(j1*a22+j2);

  // return the matrix function f(A) = f1*I + (f2-f1)*P2 [5 flops]
  a11 = diff*a11+f1;
  a22 = diff*a22+f1;
  a12 *= diff;

  // operation totals:
  // 47 flops
  // 1 sqrtf()
  // 1 fmax()
  // 1 div
  // 2 f()
  // 2 dfdx()
  
} // matrix_jacobian()

#endif // MATRIX_FUNCTION_H
