#ifndef PLASTIC_H
#define PLASTIC_H

#include <math.h>
#include <float.h> // FLT_EPSILON
#include "Elastic.h"
#include "matrix_function.h" // matrix_jacobian()

class Plastic : public Elastic {
public:

  // Plastic material parameters
  float yield; // Yield stress

  // Parameterized constructor
  Plastic(float density, float youngs_modulus, float poisson_ratio, float yield_stress) : Elastic(density, youngs_modulus, poisson_ratio) {
    yield = yield_stress;
  } // Plastic()

  // Return the number of state variables for allocation purposes
  virtual int numStateVariables(void) { return 2; }

  // Initialize the material state
#pragma acc routine seq
  virtual void initialize(float* state) {
    state[0] = 0.0f;
    state[1] = 0.0f;
  } // initialize(float* state)

  // Update the material state using the current deformation gradient F
#pragma acc routine seq
  bool update(float (&F)[2][2], float* state) {
    // compute the Lagrangian Hencky (logarithmic) strain measure: H = 0.5*log(F^T*F)
    float H11 = F[0][0]*F[0][0]+F[1][0]*F[1][0];
    float H22 = F[0][1]*F[0][1]+F[1][1]*F[1][1];
    float H12 = F[0][0]*F[0][1]+F[1][0]*F[1][1];
    float J1111,J2222,J1212,J2212,J1112,J1122;
    matrix_jacobian(H11,H22,H12,J1111,J2222,J1212,J2212,J1112,J1122,[](float lam) -> float {return 0.5f*logf(lam);},[](float lam) -> float {return 0.5f/lam;});

    // compute the dilatational and deviatoric elastic strains
    float trH = H11+H22;
    float devH11 = H11-state[0]-0.5f*trH;
    float devH12 = H12-state[1];

    // compute the effective stress
    float seff = mu2*2.0f*sqrtf(devH11*devH11+devH12*devH12);

    // compute plastic strain increment divided by Heff
    float deps = fmax(0.0f,1.0f-yield/fmax(seff,mu2*FLT_EPSILON));

    // compute the current pressure
    float press = kappa*trH;

    // update the plastic strains
    state[0] += devH11*deps;
    state[1] += devH12*deps;

    // compute the co-rotational stress tensor T = kappa*trH*I + 2*mu*dev[He-Hp]
    float mu2eff = mu2*(1.0f-deps);
    float T11 = mu2eff*devH11 + press;
    float T22 =-mu2eff*devH11 + press;
    float twoT12 = 2.0f*mu2eff*devH12;

    // transform of the co-rotational stress T into S = J : 2 T
    float S11 = 2.0f*(J1111*T11+J1122*T22+J1112*twoT12);
    float S22 = 2.0f*(J1122*T11+J2222*T22+J2212*twoT12);
    float S12 = 2.0f*(J1112*T11+J2212*T22+J1212*twoT12);

    // transform the second P-K stress into the first P-K stress: P = F*S
    // (return as F)
    float P[2][2] = { {F[0][0]*S11+F[0][1]*S12, F[0][0]*S12+F[0][1]*S22},
                      {F[1][0]*S11+F[1][1]*S12, F[1][0]*S12+F[1][1]*S22} };
    F[0][0] = P[0][0];
    F[0][1] = P[0][1];
    F[1][0] = P[1][0];
    F[1][1] = P[1][1];

    return true;
  } // update(float (&F)[2][2], float* state)

    // Update the material state using the current deformation gradient F
#pragma acc routine seq
  bool updatePlastic(float (&F)[2][2], float* state) {
    // compute the Lagrangian Hencky (logarithmic) strain measure: H = 0.5*log(F^T*F)
    float H11 = F[0][0]*F[0][0]+F[1][0]*F[1][0];
    float H22 = F[0][1]*F[0][1]+F[1][1]*F[1][1];
    float H12 = F[0][0]*F[0][1]+F[1][0]*F[1][1];
    float J1111,J2222,J1212,J2212,J1112,J1122;
    matrix_jacobian(H11,H22,H12,J1111,J2222,J1212,J2212,J1112,J1122,[](float lam) -> float {return 0.5f*logf(lam);},[](float lam) -> float {return 0.5f/lam;});

    // compute the dilatational and deviatoric elastic strains
    float trH = H11+H22;
    float devH11 = H11-state[0]-0.5f*trH;
    float devH12 = H12-state[1];

    // compute the effective stress
    float seff = mu2*2.0f*sqrtf(devH11*devH11+devH12*devH12);

    // compute plastic strain increment divided by Heff
    float deps = fmax(0.0f,1.0f-yield/fmax(seff,mu2*FLT_EPSILON));

    // compute the current pressure
    float press = kappa*trH;

    // update the plastic strains
    state[0] += devH11*deps;
    state[1] += devH12*deps;

    // compute the co-rotational stress tensor T = kappa*trH*I + 2*mu*dev[He-Hp]
    float mu2eff = mu2*(1.0f-deps);
    float T11 = mu2eff*devH11 + press;
    float T22 =-mu2eff*devH11 + press;
    float twoT12 = 2.0f*mu2eff*devH12;

    // transform of the co-rotational stress T into S = J : T
    float S11 = J1111*T11+J1122*T22+J1112*twoT12;
    float S22 = J1122*T11+J2222*T22+J2212*twoT12;
    float S12 = J1112*T11+J2212*T22+J1212*twoT12;

    // transform the second P-K stress into the first P-K stress: P = F*S
    // (return as F)
    float P[2][2] = { {F[0][0]*S11+F[0][1]*S12, F[0][0]*S12+F[0][1]*S22},
                      {F[1][0]*S11+F[1][1]*S12, F[1][0]*S12+F[1][1]*S22} };
    F[0][0] = P[0][0];
    F[0][1] = P[0][1];
    F[1][0] = P[1][0];
    F[1][1] = P[1][1];

    return true;
  } // update(float (&F)[2][2], float* state)
  
}; // Plastic()

#endif // PLASTIC_H
