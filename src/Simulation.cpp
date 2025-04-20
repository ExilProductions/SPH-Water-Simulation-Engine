#include "WaterSim/Simulation.h"
#include <cmath>
#include <iostream>
#include <random>

namespace WaterSim {

Simulation::Simulation() :
    solver_(std::make_unique<Solver>()),
    initialized_(false)
{
}

Simulation::~Simulation() {
}

//initialize solver and set params
void Simulation::initialize(const SimulationParameters& params) {
    solver_->setParameters(params);
    solver_->initialize();
    initialized_ = true;
}

//run simulation steps
void Simulation::run(int steps) {
    if (!initialized_) { //if not initialized
        std::cerr << "Simulation must be initialized before running." << std::endl;
        return;
    }
    
    for (int i = 0; i < steps; ++i) {
        step(); //call step
    }
}

//single simulation step
void Simulation::step() {
    if (!initialized_) { //if not initialized
        std::cerr << "Simulation must be initialized before stepping." << std::endl;
        return;
    }
    
    solver_->step(); //call solver step
    
    if (renderer_) { //if renderer set
        renderer_->render(solver_->getParticles()); //render particles
    }
}

//set the renderer
void Simulation::setRenderer(std::shared_ptr<IRenderer> renderer) {
    renderer_ = renderer;
    
    if (renderer_ && initialized_) { //if renderer and initialized
        const auto& params = solver_->getParameters();
        int width = static_cast<int>((params.x_max - params.x_min) * 800); //calculate width
        int height = static_cast<int>((params.y_max - params.y_min) * 800); //calculate height
        renderer_->initialize(width, height); //initialize renderer
    }
}

//get solver reference
Solver& Simulation::getSolver() {
    return *solver_;
}

//update simulation, same as step
void Simulation::update() {
    step();
}


//add boundary particles for 2D
void Simulation::addBoundaryParticles(float spacing) {
    const auto& params = solver_->getParameters();
    
    std::vector<Particle> boundary_particles;
    
    //add particles to all boundary sides
    for (float x = params.x_min; x <= params.x_max; x += spacing) {
        for (int layer = 0; layer < 3; ++layer) {
            float y = params.y_min - layer * spacing;
            Particle p(x, y);
            p.is_boundary = true;
            boundary_particles.push_back(p);
        }
    }
    
    for (float y = params.y_min; y <= params.y_max; y += spacing) {
        for (int layer = 0; layer < 3; ++layer) {
            float x = params.x_min - layer * spacing;
            Particle p(x, y);
            p.is_boundary = true;
            boundary_particles.push_back(p);
        }
    }
    
    for (float y = params.y_min; y <= params.y_max; y += spacing) {
        for (int layer = 0; layer < 3; ++layer) {
            float x = params.x_max + layer * spacing;
            Particle p(x, y);
            p.is_boundary = true;
            boundary_particles.push_back(p);
        }
    }
    
    solver_->addParticles(boundary_particles); //add boundary particles
}

//add boundary particles for 3D
void Simulation::add3DBoundaryParticles(float spacing) {
    const auto& params = solver_->getParameters();
    
    std::vector<Particle> boundary_particles;
    
    //add boundary particles for each 3D side
    for (float x = params.x_min; x <= params.x_max; x += spacing) {
        for (float z = params.z_min; z <= params.z_max; z += spacing) {
            for (int layer = 0; layer < 2; ++layer) {
                float y = params.y_min - layer * spacing;
                Particle p(x, y, z, 1.0f);
                p.is_boundary = true;
                boundary_particles.push_back(p);
            }
        }
    }
    
    for (float y = params.y_min; y <= params.y_max; y += spacing) {
        for (float z = params.z_min; z <= params.z_max; z += spacing) {
            for (int layer = 0; layer < 2; ++layer) {
                float x = params.x_min - layer * spacing;
                Particle p(x, y, z, 1.0f);
                p.is_boundary = true;
                boundary_particles.push_back(p);
            }
        }
    }
    
    for (float y = params.y_min; y <= params.y_max; y += spacing) {
        for (float z = params.z_min; z <= params.z_max; z += spacing) {
            for (int layer = 0; layer < 2; ++layer) {
                float x = params.x_max + layer * spacing;
                Particle p(x, y, z, 1.0f);
                p.is_boundary = true;
                boundary_particles.push_back(p);
            }
        }
    }
    
    for (float x = params.x_min; x <= params.x_max; x += spacing) {
        for (float y = params.y_min; y <= params.y_max; y += spacing) {
            for (int layer = 0; layer < 2; ++layer) {
                float z = params.z_min - layer * spacing;
                Particle p(x, y, z, 1.0f);
                p.is_boundary = true;
                boundary_particles.push_back(p);
            }
        }
    }
    
    for (float x = params.x_min; x <= params.x_max; x += spacing) {
        for (float y = params.y_min; y <= params.y_max; y += spacing) {
            for (int layer = 0; layer < 2; ++layer) {
                float z = params.z_max + layer * spacing;
                Particle p(x, y, z, 1.0f);
                p.is_boundary = true;
                boundary_particles.push_back(p);
            }
        }
    }
    
    solver_->addParticles(boundary_particles); //add 3D boundary particles
}

}
