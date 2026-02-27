#include "spatialHashing.h"

#include "quicksort.h"
#include "sectors.h"
#include "config.h"
#include "particle.h"

#include <string.h>
#include <math.h>
#include <stdio.h>


static int cellX, cellY;
static int i,j;

int getCellKey(int cellX, int cellY, SimConfig* config) {
    // return (cellX * config->primex + cellY * config->primey) % config->numCells;
    return cellX * config->numCellsY + cellY; // Simple row-major ordering
}

void updateCellKey(Sector* sectors, Particle* particle, SimConfig* config) {
    // Calculate cell indices based on particle position and cell size
    cellX = (int)((particle->x + config->boxLength / 2) / config->cellSize);
    cellY = (int)((particle->y + config->boxLength / 2) / config->cellSize);

    // Get list of cellKeys for neighbouring cells
    for (i = -1; i <= 1; i++) {
        for (j = -1; j <= 1; j++) {
            int neighbourX = cellX + i;
            int neighbourY = cellY + j;
            if (neighbourX >= 0 && neighbourX < config->numCellsX && neighbourY >= 0 && neighbourY < config->numCellsY) { // Catch boundary conditions
                particle->cellSearchKey[(i+1)*3 + (j+1)] = getCellKey(neighbourX, neighbourY, config);
            } else {
                particle->cellSearchKey[(i+1)*3 + (j+1)] = -1;
            }
        }
    }
    
    // Update particle's spatial hashing data
    particle->cellX = cellX;
    particle->cellY = cellY;
    particle->cellKey = getCellKey(cellX, cellY, config);

    // Update sector population count for the particle's new cell
    int sectorIndex = particle->cellKey;
    if (sectors[sectorIndex].population == 0) {
        sectors[sectorIndex].key = particle->cellKey; // Set the sector key if it's the first particle in that sector
    }
    sectors[sectorIndex].population++;
}

// void updateSpatialLookup(Sector* sectors, Particle particles[], SimConfig* config) {
//     // Reset startIndices
//     for (int i = 0; i < config->numCells; i++) {
//         sectors[i].startIndex = -1;
//     }

//     // Sort the particles by cellKey
//     quickSort(particles, 0, config->numParticles - 1);
    
//     // Get the start index for each unique cellKey
//     for (i = 0; i < config->numParticles; i++) {

//         // If the current key[i] is different from the previous key, i is the start index for that key[i]
//         // Special case for the first particle, which is always the start index for the first key.
//         if (i == 0 || particles[i - 1].cellKey != particles[i].cellKey) {
//             sectors[particles[i].cellKey].startIndex = i;
//         }
//     }
// }

void updateSpatialLookup(Sector* sectors, Particle particles[], SimConfig* config) {
    // Reset startIndices
    for (int i = 0; i < config->numCells; i++) {
        sectors[i].startIndex = -1;
    }

    // Update sector start indices
    for (int i=0; i < config->numCells; i++) {
        if (i==0) {
            sectors[i].startIndex = 0;
        } else {
            sectors[i].startIndex = sectors[i-1].startIndex + sectors[i-1].population;
            sectors[i].count = 0; // Reset counts
        }
    }

    // Insert Particles into Sorted List
    Particle sortedParticles[config->numParticles];
    for (int p=0; p < config->numParticles; p++) {
        int sectorIndex = particles[p].cellKey;
        int insertIndex = sectors[sectorIndex].startIndex + sectors[sectorIndex].count;
        sortedParticles[insertIndex] = particles[p];
        sectors[sectorIndex].count++;
    }

    // Copy sorted particles back to original array
    memcpy(particles, sortedParticles, config->numParticles * sizeof(Particle));
}