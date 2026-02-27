#include "particle.h"

Particle createParticle(){
    Particle p;
    p.ID = 0;
    p.mass = 1.0;
    p.sigma = 1.0;
    p.charge = 1.0;

    // Motion
    p.x = 0.0;
    p.y = 0.0;
    p.z = 1.0;
    p.vx = 0.0;
    p.vy = 0.0;
    p.vz = 0.0;
    p.vhx = 0.0;
    p.vhy = 0.0;
    p.vhz = 0.0;
    p.ax = 0.0;
    p.ay = 0.0;
    p.az = 0.0;
    p.Fx = 0.0;
    p.Fy = 0.0;
    p.Fz = 0.0;

    // Spatial Hashing
    p.cellKey = -1;
    p.cellX = -1;
    p.cellY = -1;
    for (int i = 0; i < 9; i++) {
        p.cellSearchKey[i] = -1;
    }

    p.isPlaced = false;
    p.isTouching = false;
    p.isEscaped = false;

    return p;
}