#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "simulation.h" // Simulation functions
#include "placement.h"  // Particle placement algorithm
#include "spatialHashing.h" // Spatial hashing functions
#include "forces.h"    // Force calculation functions

#include "config.h"     // Configuration settings
#include "particle.h"   // Particle structure

//// Constants ////
#define PI 3.14159265359f


bool DEBUG_MODE = false; // Global debug mode flag

int main(int argc, char* argv[]) {
    float clock_start, clock_end;
    clock_start = clock();

    //// System Parameters
    printf("Setting up configuration...\n");
    SimConfig config = defaultConfig();
    if (argc > 1 && strcmp(argv[1], "--debug") == 0) {
        DEBUG_MODE = true;
    }

    //// Particle Parameters
    printf("Creating particles...\n");
    Particle particles[config.numParticles];
    for (int i = 0; i < config.numParticles; i++) {
        // Generate default particle
        particles[i] = createParticle();
        // Set unique ID and particle properties
        particles[i].ID = i;
        particles[i].mass = 1.0;
        particles[i].sigma = config.sigma;
        particles[i].charge = 1.0;
    }


    //// Simulation
    // srand(42); // Seed for reproducibility
    srand(time(NULL)); // Seed by time for actual randomness
    printf("Placing particles...\n");

    placeParticles(particles, &config);


    for (int i = 0; i < config.numParticles; i++) {
        applyReflectingBoundary(&particles[i], &config); // Ensure no particles are placed outside the boundaries
        particles[i].cellKey = getCellKey((int)((particles[i].x + config.boxLength / 2) / config.cellSize), (int)((particles[i].y + config.boxLength / 2) / config.cellSize), &config); // Set initial cell key for each particle
    }
    printf("Running simulation...\n");
    runSimulation(particles, &config);

    printf("\nSimulation complete!\n");

    clock_end = clock();
    double elapsed_time = (clock_end - clock_start) / CLOCKS_PER_SEC;
    printf("Elapsed Time: %f seconds\n", elapsed_time);


    
    return 0;
}