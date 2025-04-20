#pragma once

#include "Particle.h"
#include <vector>
#include <unordered_map>
#include <functional>

namespace WaterSim {

/**
*@class SpatialGrid
*@brief grid for finding neighbors faster in 2D and 3D
*/
class SpatialGrid {
public:
    SpatialGrid(float cell_size);
    
    //add particle to grid
    void addParticle(int particle_idx, const Particle& particle);
    
    //clear grid (called before each update)
    void clear();
    
    //get nearby particles within radius
    std::vector<int> getNeighbors(const Particle& particle, float radius) const;
    
    //find particle pairs within range
    void findPairs(const std::vector<Particle>& particles, float radius, 
                  std::function<void(int, int)> callback) const;
    
private:
    float cell_size_;
    
    //map cell to list of particles
    std::unordered_map<int, std::vector<int>> grid_;
    
    //hash function for 3D cell index
    int hashFunction(int cell_x, int cell_y, int cell_z) const;
    
    //get grid cell for position
    void getCellCoords(float x, float y, float z, int& cell_x, int& cell_y, int& cell_z) const;
};

}
