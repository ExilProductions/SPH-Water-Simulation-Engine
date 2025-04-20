#include "WaterSim/SpatialGrid.h"
#include <cmath>
#include <algorithm>

namespace WaterSim {

SpatialGrid::SpatialGrid(float cell_size) : cell_size_(cell_size) {
}

void SpatialGrid::addParticle(int particle_idx, const Particle& particle) {
    int cell_x, cell_y, cell_z;
    getCellCoords(particle.x, particle.y, particle.z, cell_x, cell_y, cell_z);
    int cell_hash = hashFunction(cell_x, cell_y, cell_z);
    grid_[cell_hash].push_back(particle_idx);
}

void SpatialGrid::clear() {
    grid_.clear();
}

std::vector<int> SpatialGrid::getNeighbors(const Particle& particle, float radius) const {
    std::vector<int> neighbors;
    
    int cell_radius = static_cast<int>(std::ceil(radius / cell_size_));
    int center_cell_x, center_cell_y, center_cell_z;
    getCellCoords(particle.x, particle.y, particle.z, center_cell_x, center_cell_y, center_cell_z);
    
    
    for (int i = -cell_radius; i <= cell_radius; ++i) {
        for (int j = -cell_radius; j <= cell_radius; ++j) {
            for (int k = -cell_radius; k <= cell_radius; ++k) {
                int cell_x = center_cell_x + i;
                int cell_y = center_cell_y + j;
                int cell_z = center_cell_z + k;
                int cell_hash = hashFunction(cell_x, cell_y, cell_z);
                
                auto it = grid_.find(cell_hash);
                if (it == grid_.end())
                    continue;
                    
                for (int particle_idx : it->second) {
                    neighbors.push_back(particle_idx);
                }
            }
        }
    }
    
    return neighbors;
}

void SpatialGrid::findPairs(const std::vector<Particle>& particles, float radius, 
                           std::function<void(int, int)> callback) const {
    float radius_squared = radius * radius;
    
    
    for (const auto& cell_entry : grid_) {
        const auto& particles_in_cell = cell_entry.second;
        
        
        for (size_t i = 0; i < particles_in_cell.size(); ++i) {
            int idx_i = particles_in_cell[i];
            
            
            for (size_t j = i + 1; j < particles_in_cell.size(); ++j) {
                int idx_j = particles_in_cell[j];
                
                if (particles[idx_i].distanceSquaredTo(particles[idx_j]) <= radius_squared) {
                    callback(idx_i, idx_j);
                }
            }
            
            
            int cell_x, cell_y, cell_z;
            getCellCoords(particles[idx_i].x, particles[idx_i].y, particles[idx_i].z, 
                         cell_x, cell_y, cell_z);
            
            
            
            const int neighbor_cells[13][3] = {
                {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {-1, 1, 0},
                {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}, {-1, 1, 1}, {-1, 0, 1}, {-1, -1, 1}, {0, -1, 1}, {1, -1, 1}
            };
            
            for (int n = 0; n < 13; ++n) {
                int neighbor_x = cell_x + neighbor_cells[n][0];
                int neighbor_y = cell_y + neighbor_cells[n][1];
                int neighbor_z = cell_z + neighbor_cells[n][2];
                int neighbor_hash = hashFunction(neighbor_x, neighbor_y, neighbor_z);
                
                auto it = grid_.find(neighbor_hash);
                if (it == grid_.end())
                    continue;
                    
                for (int idx_j : it->second) {
                    if (particles[idx_i].distanceSquaredTo(particles[idx_j]) <= radius_squared) {
                        callback(idx_i, idx_j);
                    }
                }
            }
        }
    }
}

int SpatialGrid::hashFunction(int cell_x, int cell_y, int cell_z) const {
    
    
    return 541 * cell_x + 79 * cell_y + 31 * cell_z;
}

void SpatialGrid::getCellCoords(float x, float y, float z, int& cell_x, int& cell_y, int& cell_z) const {
    cell_x = static_cast<int>(std::floor(x / cell_size_));
    cell_y = static_cast<int>(std::floor(y / cell_size_));
    cell_z = static_cast<int>(std::floor(z / cell_size_));
}

} 