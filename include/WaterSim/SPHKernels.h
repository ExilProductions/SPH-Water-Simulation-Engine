#pragma once

#include <cmath>

namespace WaterSim {

/**
*@class SPHKernels
*@brief Collection of smoothing kernels for SPH calculations in 2D and 3D
*/
class SPHKernels {
public:
    //init with smoothing radius
    SPHKernels(float h);
    
    //poly6 kernel for density
    float poly6(float r_squared) const;
    
    //gradient of poly6
    float poly6Gradient(float r, float r_squared) const;
    
    //spiky kernel for pressure force
    float spiky(float r) const;
    
    //gradient of spiky
    float spikyGradient(float r) const;
    
    //viscosity kernel
    float viscosity(float r) const;
    
    //viscosity laplacian
    float viscosityLaplacian(float r) const;
    
    //get smoothing radius
    float getSmoothingRadius() const { return h_; }
    
    //get squared smoothing radius
    float getSmoothingRadiusSquared() const { return h_squared_; }
    
    //set 3D kernels usage
    void setUse3D(bool use_3d) { use_3d_ = use_3d; }
    
    //check 3D kernels
    bool isUsing3D() const { return use_3d_; }
    
private:
    float h_;         //smoothing radius
    float h_squared_; //squared smoothing radius
    float h_cubed_;   //cubed smoothing radius
    float h_fourth_;  //h^4
    float h_fifth_;   //h^5
    float h_sixth_;   //h^6
    float h_eighth_;  //h^8
    float h_ninth_;   //h^9
    
    bool use_3d_;     //whether to use 3D kernels
    
    //precomputed coefficients
    float poly6_coeff_2d_;
    float poly6_coeff_3d_;
    float spiky_coeff_2d_;
    float spiky_coeff_3d_;
    float viscosity_coeff_2d_;
    float viscosity_coeff_3d_;
};

} // namespace WaterSim
