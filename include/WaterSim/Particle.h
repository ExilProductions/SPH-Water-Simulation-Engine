#pragma once

#include <vector>
#include <cmath>

namespace WaterSim {

/**
*@class Particle
*@brief Represents a single particle in the SPH simulation
*/
class Particle {
public:
    //core properties
    float x, y, z;           //position (3D)
    float vx, vy, vz;        //velocity (3D)
    float ax, ay, az;        //acceleration (3D)
    float density;           //fluid density at this particle
    float pressure;          //pressure at this particle
    float mass;              //particle mass
    
    //extra properties for sim
    float predicted_density;
    int id;                  //unique particle id
    bool is_boundary;        //flag for boundary particles
    
    Particle() : 
        x(0.0f), y(0.0f), z(0.0f),
        vx(0.0f), vy(0.0f), vz(0.0f),
        ax(0.0f), ay(0.0f), az(0.0f),
        density(0.0f), pressure(0.0f),
        mass(1.0f), predicted_density(0.0f), 
        id(0), is_boundary(false) {}
    
    Particle(float x, float y, float mass = 1.0f) : 
        x(x), y(y), z(0.0f),
        vx(0.0f), vy(0.0f), vz(0.0f),
        ax(0.0f), ay(0.0f), az(0.0f),
        density(0.0f), pressure(0.0f),
        mass(mass), predicted_density(0.0f), 
        id(0), is_boundary(false) {}
    
    Particle(float x, float y, float z, float mass = 1.0f) : 
        x(x), y(y), z(z),
        vx(0.0f), vy(0.0f), vz(0.0f),
        ax(0.0f), ay(0.0f), az(0.0f),
        density(0.0f), pressure(0.0f),
        mass(mass), predicted_density(0.0f), 
        id(0), is_boundary(false) {}
    
    //calc distance to another particle
    float distanceTo(const Particle& other) const {
        float dx = x - other.x;
        float dy = y - other.y;
        float dz = z - other.z;
        return std::sqrt(dx*dx + dy*dy + dz*dz);
    }
    
    //calc squared distance to another particle
    float distanceSquaredTo(const Particle& other) const {
        float dx = x - other.x;
        float dy = y - other.y;
        float dz = z - other.z;
        return dx*dx + dy*dy + dz*dz;
    }
};

}
