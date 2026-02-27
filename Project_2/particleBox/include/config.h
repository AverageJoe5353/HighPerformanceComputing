#ifndef CONFIG_H
#define CONFIG_H

#include <stdbool.h>

typedef struct {
    // Simulation Parameters
    int numParticles;    // Number of particles to simulate
    double boxLength;    // Side length of box
    double sigma;        // Particle diameter
    double epsilon;      // WCA Adjustment
    double kT;           // System temperature [kT]
    double gamma;        // Friction Coefficient
    double timeStep;     // Time Step Size
    double maxStep;      // Max number of time steps
    double r_cutoff;     // Interaction cutoff distance
    int recordInterval; // Interval for recording trajectory data
    bool recordToggle; // Toggle for recording trajectory data

    double V; // Voltage across nanopore
    double pore_a; // Half-width of nanopore

    int numCells;     // Number of cells in spatial hashing grid
    double cellSize;  // Size of each cell in spatial hashing grid
    int numCellsX;    // Number of cells in X direction
    int numCellsY;    // Number of cells in Y direction
    int primex;       // Prime number for spatial hashing key generation
    int primey;       // Prime number for spatial hashing key generation

} SimConfig;

SimConfig defaultConfig(); // Default simulation configuration
void debugrep(const char* message); // Debug print function
void debugstage(const char* message); // Debug stage function

#endif // CONFIG_H