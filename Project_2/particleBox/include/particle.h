#ifndef PARTICLE_H
#define PARTICLE_H

#include <stdbool.h>

// Particle structure definition
typedef struct {
    // Particle Attributes
    int ID;    // Particle ID
    double mass;  // Particle Mass
    double sigma; // Particle Size
    double charge; // Particle Charge

    // Motion
    double x;      // x position
    double y;      // y position
    double z;      // z position
    double vx;     // x velocity
    double vy;     // y velocity
    double vz;     // z velocity
    double vhx;    // x half-step velocity
    double vhy;    // y half-step velocity
    double vhz;    // z half-step velocity
    double ax;     // x acceleration
    double ay;     // y acceleration
    double az;     // z acceleration
    double Fx;     // x net force
    double Fy;     // y net force
    double Fz;     // z net force

    // Spatial Hashing
    int cellKey; // Key for spatial hashing grid cell
    int cellX;   // Cell X index
    int cellY;   // Cell Y index
    int cellSearchKey[9]; // Indices of cells for neighbour search

    // Status Flagging
    bool isPlaced;   // Flag for placement status
    bool isTouching; // Flag for overlap status
    bool isEscaped;  // Flag for breaking boundaries
} Particle;

Particle createParticle(); // Generates blank prototype particle struct


#endif // PARTICLE_H