#ifndef PLASTIC_H
#define PLASTIC_H

#include <math.h>
#include "Elastic.h"
#include "matrix_function.h"

class Plastic : Elastic {

  float yield; // Yield stress

  Plastic(float density, float youngs_modulus, float poisson_ratio, float yield_stress) : Elastic(float density, float youngs_modulus, float poisson_ratio) {
    yield = yield_stress;
  } // Plastic()

  virtual bool update(float** F, float* state) {
    // compute the Eulerian Hencky (logarithmic) strain measure: H = 0.5*log(F*F^T)
    float H11 = F[0][0]*F[0][0]+F[0][1]*F[0][1];
    float H22 = F[1][0]*F[1][0]+F[1][1]*F[1][1];
    float H12 = F[0][0]*F[1][0]+F[0][1]*F[1][1];
    matrix_function(H11,H22,H12,[](float lam) -> float {return 0.5f*logf(lam);});

    // compute the dilatational and deviatoric elastic strains
    float trH = H11+H22;
    float devH11 = H11-state_vars[n_state_vars*(4*ie+q)  ]-0.5f*trH;
    float devH12 = H12-state_vars[n_state_vars*(4*ie+q)+1];

    // compute the effective stress
    float seff = mu2*2.0f*sqrtf(devH11*devH11+devH12*devH12);

    // compute plastic strain increment divided by Heff
    float deps = fmax(0.0f,1.0f-yield/(seff+mu2*1.0e-16));

    // compute the current pressure
    float press = kappa*trH;

    // update the plastic strains
    state[0] += devH11*deps;
    state[1] += devH12*deps;

    // compute the Kirchhoff stress tensor
    float mu2eff = mu2*(1.0f-deps);
    float tau11 = mu2eff*devH11 + press;
    float tau22 =-mu2eff*devH11 + press;
    float tau12 = mu2eff*devH12;

    // compute the Kirchhoff stress tensor
    // for a Hencky elastic model: tau = 2*mu*H + lam*tr(H)*I
    //float tau11 = pmod*H11 + lam*H22;
    //float tau22 = pmod*H22 + lam*H11;
    //float tau12 = mu2*H12;

    // transform the Kirchhoff stress into the first P-K stress: P = tau*F^-T
    // (times the differential volume dV)
    float dVdetF = dV/(F[0][0]*F[1][1] - F[0][1]*F[1][0]);
    float P_dV[2][2] = { {(+tau11*F[1][1]-tau12*F[0][1])*dVdetF,
			  (-tau11*F[1][0]+tau12*F[0][0])*dVdetF},
			 {(+tau12*F[1][1]-tau22*F[0][1])*dVdetF,
			  (-tau12*F[1][0]+tau22*F[0][0])*dVdetF} };

    // compute the Kirchhoff stress tensor
    // for a Hencky elastic model: tau = 2*mu*H + lam*tr(H)*I
    float tau11 = pmod*H11 + lam*H22;
    float tau22 = pmod*H22 + lam*H11;
    float tau12 = mu2*H12;

    // transform the Kirchhoff stress into the first P-K stress: P = tau*F^-T
    // (return as F)
    float idetF = 1.0/(F[0][0]*F[1][1] - F[0][1]*F[1][0]);
    F = { {(+tau11*F[1][1]-tau12*F[0][1])*idetF,
           (-tau11*F[1][0]+tau12*F[0][0])*idetF},
          {(+tau12*F[1][1]-tau22*F[0][1])*idetF,
           (-tau12*F[1][0]+tau22*F[0][0])*idetF} };

    return true;
  } // update(float** F, float* state)
  
} // Plastic()

#endif // PLASTIC_H
