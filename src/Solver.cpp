#include "WaterSim/Solver.h"
#include <algorithm>
#include <iostream>

namespace WaterSim {

Solver::Solver(const SimulationParameters& params) :
    grid_(params.smoothing_radius),
    kernels_(params.smoothing_radius),
    params_(params)
{
    kernels_.setUse3D(params.use_3d);
}

void Solver::initialize() {
    
    grid_ = SpatialGrid(params_.smoothing_radius);
    
    
    kernels_ = SPHKernels(params_.smoothing_radius);
    kernels_.setUse3D(params_.use_3d);
    
    
    particles_.clear();
}

int Solver::addParticle(float x, float y, float z, bool is_boundary) {
    Particle p(x, y, 0.0f, params_.particle_mass); 
    p.z = z;
    p.id = static_cast<int>(particles_.size());
    p.is_boundary = is_boundary;
    particles_.push_back(p);
    return p.id;
}

void Solver::addParticles(const std::vector<Particle>& new_particles) {
    size_t start_idx = particles_.size();
    particles_.insert(particles_.end(), new_particles.begin(), new_particles.end());
    
    
    for (size_t i = 0; i < new_particles.size(); ++i) {
        particles_[start_idx + i].id = static_cast<int>(start_idx + i);
    }
}

void Solver::removeParticle(int particle_id) {
    if (particle_id >= 0 && particle_id < static_cast<int>(particles_.size())) {
        
        std::swap(particles_[particle_id], particles_.back());
        particles_.pop_back();
        
        
        if (particle_id < static_cast<int>(particles_.size())) {
            particles_[particle_id].id = particle_id;
        }
    }
}

void Solver::clear() {
    particles_.clear();
}

void Solver::step() {
    updateGrid();
    computeDensityPressure();
    predictDensities();
    solveIncompressibility();
    computeForces();
    integrate();
    updateBoundaries();
}

void Solver::setParameters(const SimulationParameters& params) {
    params_ = params;
    kernels_ = SPHKernels(params_.smoothing_radius);
    kernels_.setUse3D(params_.use_3d);
    grid_ = SpatialGrid(params_.smoothing_radius);
}

const SimulationParameters& Solver::getParameters() const {
    return params_;
}

const std::vector<Particle>& Solver::getParticles() const {
    return particles_;
}

void Solver::updateGrid() {
    grid_.clear();
    for (size_t i = 0; i < particles_.size(); ++i) {
        grid_.addParticle(i, particles_[i]);
    }
}

void Solver::computeDensityPressure() {
    float h_squared = kernels_.getSmoothingRadiusSquared();
    
    
    for (size_t i = 0; i < particles_.size(); ++i) {
        Particle& pi = particles_[i];
        
        
        pi.density = kernels_.poly6(0.0f) * pi.mass;
        
        
        std::vector<int> neighbors = grid_.getNeighbors(pi, params_.smoothing_radius);
        
        
        for (int j : neighbors) {
            if (static_cast<int>(i) == j) continue;
            
            Particle& pj = particles_[j];
            float r_squared = pi.distanceSquaredTo(pj);
            
            if (r_squared < h_squared) {
                pi.density += pj.mass * kernels_.poly6(r_squared);
            }
        }
        
        
        pi.pressure = params_.gas_constant * (pi.density - params_.rest_density);
        if (pi.pressure < 0.0f) pi.pressure = 0.0f; 
    }
}

void Solver::predictDensities() {
    
    for (size_t i = 0; i < particles_.size(); ++i) {
        Particle& pi = particles_[i];
        
        
        pi.predicted_density = pi.density;
        
        std::vector<int> neighbors = grid_.getNeighbors(pi, params_.smoothing_radius);
        
        for (int j : neighbors) {
            if (static_cast<int>(i) == j) continue;
            
            Particle& pj = particles_[j];
            float r = pi.distanceTo(pj);
            
            if (r < params_.smoothing_radius) {
                
                float dx = pi.x - pj.x;
                float dy = pi.y - pj.y;
                float dz = pi.z - pj.z;
                float dvx = pi.vx - pj.vx;
                float dvy = pi.vy - pj.vy;
                float dvz = pi.vz - pj.vz;
                
                
                float div = (dx * dvx + dy * dvy + dz * dvz) / (r * r + 0.0001f);
                
                
                pi.predicted_density -= div * params_.time_step * kernels_.spiky(r) * pj.mass;
            }
        }
    }
}

void Solver::solveIncompressibility() {
    
    for (int iter = 0; iter < params_.pressure_iterations; ++iter) {
        for (size_t i = 0; i < particles_.size(); ++i) {
            if (particles_[i].is_boundary) continue; 
            
            Particle& pi = particles_[i];
            std::vector<int> neighbors = grid_.getNeighbors(pi, params_.smoothing_radius);
            
            float density_constraint = pi.predicted_density - params_.rest_density;
            if (density_constraint <= 0.0f) continue;
            
            float sum_gradient_squared = 0.0f;
            std::vector<std::pair<int, float>> grad_constraints;
            
            for (int j : neighbors) {
                if (static_cast<int>(i) == j) continue;
                
                Particle& pj = particles_[j];
                float r = pi.distanceTo(pj);
                
                if (r < params_.smoothing_radius) {
                    float grad = kernels_.spikyGradient(r);
                    grad_constraints.push_back(std::make_pair(j, grad));
                    sum_gradient_squared += grad * grad;
                }
            }
            
            
            if (sum_gradient_squared < 0.0001f) continue;
            
            
            float lambda = -density_constraint / (sum_gradient_squared + params_.cfm_epsilon);
            
            
            for (const auto& grad_pair : grad_constraints) {
                int j = grad_pair.first;
                float grad = grad_pair.second;
                
                Particle& pj = particles_[j];
                
                
                if (pj.is_boundary) continue;
                
                float dx = pi.x - pj.x;
                float dy = pi.y - pj.y;
                float dz = pi.z - pj.z;
                float r = std::sqrt(dx*dx + dy*dy + dz*dz) + 0.0001f;
                
                float factor = lambda * grad / r;
                
                
                float corr_x = dx * factor;
                float corr_y = dy * factor;
                float corr_z = dz * factor;
                
                pi.x += corr_x * 0.5f;
                pi.y += corr_y * 0.5f;
                pi.z += corr_z * 0.5f;
                pj.x -= corr_x * 0.5f;
                pj.y -= corr_y * 0.5f;
                pj.z -= corr_z * 0.5f;
            }
        }
    }
}

void Solver::applyPressureForces() {
    for (size_t i = 0; i < particles_.size(); ++i) {
        Particle& pi = particles_[i];
        if (pi.is_boundary) continue; 
        
        std::vector<int> neighbors = grid_.getNeighbors(pi, params_.smoothing_radius);
        float pressure_force_x = 0.0f;
        float pressure_force_y = 0.0f;
        float pressure_force_z = 0.0f;
        
        for (int j : neighbors) {
            if (static_cast<int>(i) == j) continue;
            
            Particle& pj = particles_[j];
            float r = pi.distanceTo(pj);
            
            if (r < params_.smoothing_radius && r > 0.0001f) {
                float dx = pi.x - pj.x;
                float dy = pi.y - pj.y;
                float dz = pi.z - pj.z;
                
                
                float pressure_term = (pi.pressure + pj.pressure) / (2.0f * pj.density);
                float force_magnitude = -pj.mass * pressure_term * kernels_.spikyGradient(r);
                
                pressure_force_x += force_magnitude * dx / r;
                pressure_force_y += force_magnitude * dy / r;
                pressure_force_z += force_magnitude * dz / r;
            }
        }
        
        pi.ax += pressure_force_x;
        pi.ay += pressure_force_y;
        pi.az += pressure_force_z;
    }
}

void Solver::applyViscosityForces() {
    for (size_t i = 0; i < particles_.size(); ++i) {
        Particle& pi = particles_[i];
        if (pi.is_boundary) continue; 
        
        std::vector<int> neighbors = grid_.getNeighbors(pi, params_.smoothing_radius);
        float visc_force_x = 0.0f;
        float visc_force_y = 0.0f;
        float visc_force_z = 0.0f;
        
        for (int j : neighbors) {
            if (static_cast<int>(i) == j) continue;
            
            Particle& pj = particles_[j];
            float r = pi.distanceTo(pj);
            
            if (r < params_.smoothing_radius) {
                
                float dvx = pj.vx - pi.vx;
                float dvy = pj.vy - pi.vy;
                float dvz = pj.vz - pi.vz;
                
                
                float visc_term = params_.viscosity * pj.mass / pj.density;
                float visc_kernel = kernels_.viscosityLaplacian(r);
                
                visc_force_x += visc_term * dvx * visc_kernel;
                visc_force_y += visc_term * dvy * visc_kernel;
                visc_force_z += visc_term * dvz * visc_kernel;
            }
        }
        
        pi.ax += visc_force_x;
        pi.ay += visc_force_y;
        pi.az += visc_force_z;
    }
}

void Solver::computeForces() {
    
    for (auto& p : particles_) {
        p.ax = params_.gravity_x;
        p.ay = params_.gravity_y;
        p.az = params_.gravity_z;
    }
    
    
    applyPressureForces();
    
    
    applyViscosityForces();
}

void Solver::integrate() {
    float dt = params_.time_step;
    
    
    for (auto& p : particles_) {
        if (p.is_boundary) continue; 
        
        
        p.vx += p.ax * dt;
        p.vy += p.ay * dt;
        p.vz += p.az * dt;
        
        
        p.x += p.vx * dt;
        p.y += p.vy * dt;
        p.z += p.vz * dt;
    }
}

void Solver::updateBoundaries() {
    float damping = params_.boundary_damping;
    
    for (auto& p : particles_) {
        if (p.is_boundary) continue; 
        
        enforceBoundaryConditions(p);
    }
}

void Solver::enforceBoundaryConditions(Particle& p) {
    
    float damping = params_.boundary_damping;
    
    
    if (p.x < params_.x_min) {
        p.vx *= -damping;
        p.x = params_.x_min;
    }
    else if (p.x > params_.x_max) {
        p.vx *= -damping;
        p.x = params_.x_max;
    }
    
    
    if (p.y < params_.y_min) {
        p.vy *= -damping;
        p.y = params_.y_min;
    }
    else if (p.y > params_.y_max) {
        p.vy *= -damping;
        p.y = params_.y_max;
    }
    
    
    if (params_.use_3d) {
        if (p.z < params_.z_min) {
            p.vz *= -damping;
            p.z = params_.z_min;
        }
        else if (p.z > params_.z_max) {
            p.vz *= -damping;
            p.z = params_.z_max;
        }
    }
}

} 