#pragma once

#include "Solver.h"
#include "Renderer.h"
#include <memory>
#include <string>

namespace WaterSim {

/**
*@class Simulation
*@brief coordinates solver and renderer
*/
class Simulation {
public:
    Simulation();
    ~Simulation();
    
    //init sim
    void initialize(const SimulationParameters& params = SimulationParameters());
    
    //run sim steps
    void run(int steps = 1);
    
    //step sim
    void step();
    
    //set renderer
    void setRenderer(std::shared_ptr<IRenderer> renderer);
    
    //get solver
    Solver& getSolver();
    
    //update and render frame
    void update();
    
private:
    std::unique_ptr<Solver> solver_;
    std::shared_ptr<IRenderer> renderer_;
    bool initialized_;
    
    //add boundary particles
    void addBoundaryParticles(float spacing);
    
    //add 3D boundary particles
    void add3DBoundaryParticles(float spacing);
};

}
