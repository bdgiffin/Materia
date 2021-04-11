#ifndef MATERIAL_H
#define MATERIAL_H

class Material {

  float rho; // Mass density

  Material(float density) {
    rho = density;
  }

  virtual bool update(float** F, float* state) = 0;
  
} // Material

#endif // MATERIAL_H
