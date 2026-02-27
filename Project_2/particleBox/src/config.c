#include "config.h"

#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>

#define N 1000       // Number of Particles
#define L 50.0     // Side length of box
#define SIGMA 1.0   // Particle diameter
#define EPSILON 1.0 // WCA Adjustment
#define KT 1.0      // System temperature
#define g 1.0       // Friction Coefficient
#define DT 0.01     // Time Step Size
#define TMAX 10000.0 // Max number of time steps

    

SimConfig defaultConfig(){
    SimConfig config;
    config.numParticles = N;
    config.boxLength = L;
    config.sigma = SIGMA;
    config.epsilon = EPSILON;
    config.kT = KT;
    config.gamma = g;
    config.timeStep = DT;
    config.maxStep = TMAX;
    config.r_cutoff = pow(2.0, 1.0/6.0) * config.sigma; // Interaction cutoff distance
    config.recordInterval = 10; // Interval for recording trajectory data
    config.recordToggle = true; // Toggle for recording trajectory data

    config.V = 10.0; // Voltage across nanopore
    config.pore_a = 0.50; // Half-width of nanopore

    
    config.cellSize = 2.0;
    config.numCells = (int)(ceil(config.boxLength / config.cellSize)) * (int)(ceil(config.boxLength / config.cellSize)); // Number of cells in spatial hashing grid, based on box length and interaction radius
    config.numCellsX = (int)(ceil(config.boxLength / config.cellSize));
    config.numCellsY = (int)(ceil(config.boxLength / config.cellSize));
    config.primex = 150503;
    config.primey = 652331;

    return config;
}

extern bool DEBUG_MODE;
void debugrep(const char* message) {
    if (DEBUG_MODE) {
        printf("\n\033[0;31m[DEBUG]\033[0m %s\n", message);
    }
}

void debugstage(const char* message) {
    if (DEBUG_MODE) {
        printf("\r\033[0;31m[DEBUG]\033[0m %s", message);
    }
}