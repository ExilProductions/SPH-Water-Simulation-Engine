#include "WaterSim/SPHKernels.h"
#include <cmath>

namespace WaterSim {

SPHKernels::SPHKernels(float h) : 
    h_(h),
    h_squared_(h * h),
    h_cubed_(h_squared_ * h),
    h_fourth_(h_squared_ * h_squared_),
    h_fifth_(h_fourth_ * h),
    h_sixth_(h_cubed_ * h_cubed_),
    h_eighth_(h_fourth_ * h_fourth_),
    h_ninth_(h_cubed_ * h_cubed_ * h_cubed_),
    use_3d_(false)
{
    poly6_coeff_2d_ = 4.0f / (M_PI * h_eighth_);
    spiky_coeff_2d_ = 10.0f / (M_PI * h_fifth_);
    viscosity_coeff_2d_ = 40.0f / (M_PI * h_squared_);
    
    poly6_coeff_3d_ = 315.0f / (64.0f * M_PI * h_ninth_);
    spiky_coeff_3d_ = 15.0f / (M_PI * h_sixth_);
    viscosity_coeff_3d_ = 15.0f / (2.0f * M_PI * h_cubed_);
}

float SPHKernels::poly6(float r_squared) const {
    if (r_squared >= h_squared_)
        return 0.0f;
        
    float term = h_squared_ - r_squared;
    float term_cubed = term * term * term;
    
    if (use_3d_) {
        return poly6_coeff_3d_ * term_cubed;
    } else {
        return poly6_coeff_2d_ * term_cubed;
    }
}

float SPHKernels::poly6Gradient(float r, float r_squared) const {
    if (r_squared >= h_squared_)
        return 0.0f;
        
    float term = h_squared_ - r_squared;
    float term_squared = term * term;
    
    if (use_3d_) {
        return -6.0f * poly6_coeff_3d_ * term_squared * r;
    } else {
        return -6.0f * poly6_coeff_2d_ * term_squared * r;
    }
}

float SPHKernels::spiky(float r) const {
    if (r >= h_)
        return 0.0f;
        
    float term = h_ - r;
    float term_cubed = term * term * term;
    
    if (use_3d_) {
        return spiky_coeff_3d_ * term_cubed;
    } else {
        return spiky_coeff_2d_ * term_cubed;
    }
}

float SPHKernels::spikyGradient(float r) const {
    if (r >= h_ || r < 0.0001f)
        return 0.0f;
        
    float term = h_ - r;
    float term_squared = term * term;
    
    if (use_3d_) {
        return -3.0f * spiky_coeff_3d_ * term_squared / r;
    } else {
        return -3.0f * spiky_coeff_2d_ * term_squared / r;
    }
}

float SPHKernels::viscosity(float r) const {
    if (r >= h_)
        return 0.0f;
    
    if (use_3d_) {
        float term1 = -r * r * r / (2.0f * h_cubed_);
        float term2 = r * r / h_squared_;
        float term3 = h_ / (2.0f * r) - 1.0f;
        
        return viscosity_coeff_3d_ * (term1 + term2 + term3);
    } else {
        float term = h_ - r;
        return viscosity_coeff_2d_ * term;
    }
}

float SPHKernels::viscosityLaplacian(float r) const {
    if (r >= h_)
        return 0.0f;
    
    if (use_3d_) {
        return viscosity_coeff_3d_ * (h_ - r);
    } else {
        return viscosity_coeff_2d_ * (h_ - r);
    }
}

}