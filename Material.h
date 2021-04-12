#ifndef MATERIAL_H
#define MATERIAL_H

class Material {
public:

  // Material parameters
  float rho; // Mass density

  // Parameterized constructor
  Material(float density) {
    rho = density;
  } // Material()

  // Return the number of state variables for allocation purposes
  virtual int numStateVariables(void) = 0;

  // Initialize the material state
#pragma acc routine seq
  virtual void initialize(float* state) = 0;

  // Update the material state using the current deformation gradient F
#pragma acc routine seq
  bool update(float (&F)[2][2], float* state) { return false; }
  
}; // Material

#endif // MATERIAL_H
