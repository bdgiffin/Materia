#ifndef ELASTIC_H
#define ELASTIC_H

#include <math.h>
#include "Material.h"
#include "matrix_function.h" // matrix_function()

class Elastic : public Material {
public:

  // Elastic material parameters
  float E;     // Young's modulus
  float nu;    // Poisson's ratio
  float mu2;   // Twice the shear modulus
  float lam;   // Lame parameter
  float kappa; // Bulk modulus
  float pmod;  // P-wave modulus
  float c;     // Sound speed

  // Parameterized constructor
  Elastic(float density, float youngs_modulus, float poisson_ratio) : Material(density) {
    E = youngs_modulus;
    nu = poisson_ratio;
    mu2 = E/(1.0+nu);
    lam = mu2*nu/(1.0-2.0*nu);
    kappa = lam+0.5f*mu2;
    pmod = lam+mu2;
    c = sqrt(pmod/rho);
  } // Elastic()

  // Return the number of state variables for allocation purposes
  virtual int numStateVariables(void) { return 0; }

  // Initialize the material state
#pragma acc routine seq
  virtual void initialize(float* state) { }

  // Update the material state using the current deformation gradient F
#pragma acc routine seq
  bool update(float (&F)[2][2], float* state) {
    // compute the Eulerian Hencky (logarithmic) strain measure: H = 0.5*log(F*F^T)
    float H11 = F[0][0]*F[0][0]+F[0][1]*F[0][1];
    float H22 = F[1][0]*F[1][0]+F[1][1]*F[1][1];
    float H12 = F[0][0]*F[1][0]+F[0][1]*F[1][1];
    matrix_function(H11,H22,H12,[](float lam) -> float {return 0.5f*logf(lam);});

    // compute the Kirchhoff stress tensor
    // for a Hencky elastic model: tau = 2*mu*H + lam*tr(H)*I
    float tau11 = pmod*H11 + lam*H22;
    float tau22 = pmod*H22 + lam*H11;
    float tau12 = mu2*H12;

    // transform the Kirchhoff stress into the first P-K stress: P = tau*F^-T
    // (return as F)
    float idetF = 1.0/(F[0][0]*F[1][1] - F[0][1]*F[1][0]);
    float P[2][2] = { {(+tau11*F[1][1]-tau12*F[0][1])*idetF,
                       (-tau11*F[1][0]+tau12*F[0][0])*idetF},
                      {(+tau12*F[1][1]-tau22*F[0][1])*idetF,
                       (-tau12*F[1][0]+tau22*F[0][0])*idetF} };
    F[0][0] = P[0][0];
    F[0][1] = P[0][1];
    F[1][0] = P[1][0];
    F[1][1] = P[1][1];
    
    return true;
  } // update(float (&F)[2][2], float* state)
  
}; // Elastic

#endif // ELASTIC_H
