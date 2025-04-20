#pragma once

#include "Particle.h"
#include "SpatialGrid.h"
#include "SPHKernels.h"
#include <vector>

namespace WaterSim {

/**
*@struct SimulationParameters
*@brief parameters to set up the simulation
*/
struct SimulationParameters {
    float time_step = 0.005f;           //time step
    float smoothing_radius = 0.1f;      //smoothing radius
    float rest_density = 1000.0f;       //rest density
    float gas_constant = 2000.0f;       //gas constant
    float viscosity = 0.1f;             //viscosity
    float gravity_x = 0.0f;             //gravity x
    float gravity_y = -9.8f;            //gravity y
    float gravity_z = 0.0f;             //gravity z
    float boundary_damping = 0.5f;      //boundary damping
    float particle_mass = 1.0f;         //particle mass
    
    //boundary settings
    float x_min = 0.0f, x_max = 1.0f;   //x bounds
    float y_min = 0.0f, y_max = 1.0f;   //y bounds
    float z_min = 0.0f, z_max = 1.0f;   //z bounds
    
    //solver stability
    float cfm_epsilon = 0.1f;           //cfm param
    int pressure_iterations = 3;        //pressure solver steps
    
    //dimensional control
    bool use_3d = false;                //3D sim flag
};

/**
*@class Solver
*@brief main SPH solver that handles the simulation
*/
class Solver {
public:
    Solver(const SimulationParameters& params = SimulationParameters());
    
    //init sim
    void initialize();
    
    //add particle to sim
    int addParticle(float x, float y, float z = 0.0f, bool is_boundary = false);
    
    //add multiple particles
    void addParticles(const std::vector<Particle>& particles);
    
    //remove particle from sim
    void removeParticle(int particle_id);
    
    //clear all particles
    void clear();
    
    //run one sim step
    void step();
    
    //set params
    void setParameters(const SimulationParameters& params);
    
    //get current params
    const SimulationParameters& getParameters() const;
    
    //get read-only particle list
    const std::vector<Particle>& getParticles() const;
    
    //update boundaries (called by step())
    void updateBoundaries();
    
private:
    std::vector<Particle> particles_;
    SpatialGrid grid_;
    SPHKernels kernels_;
    SimulationParameters params_;
    
    //main sim steps
    void updateGrid();
    void computeDensityPressure();
    void computeForces();
    void integrate();
    
    //helper methods
    void predictDensities();
    void solveIncompressibility();
    void applyPressureForces();
    void applyViscosityForces();
    
    //handle boundary conditions
    void enforceBoundaryConditions(Particle& particle);
};

}
