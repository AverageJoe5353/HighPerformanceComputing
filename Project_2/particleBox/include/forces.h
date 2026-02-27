#ifndef FORCES_H
#define FORCES_H

#include "config.h"
#include "particle.h"

void applyForce_LJ(Particle* p1, Particle* p2, SimConfig* config); // Applies Lennard-Jones force between two particles and updates their net forces
void applyReflectingBoundary(Particle* p, SimConfig* config); // Applies reflecting boundary conditions to a particle based on its position and the box length
double d1(float x, float y, float z, float a); // Distance to left edge of pore located at (-a,0)
double d2(float x, float y, float z, float a); // Distance to right edge of pore located at (a,0)
void applyForce_nanoporeEfield(Particle* p, SimConfig* config); // Applies electric field force on a particle based on its position and the nanopore configuration

#endif // FORCES_H