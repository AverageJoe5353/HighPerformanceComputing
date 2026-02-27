#ifndef SPATIAL_HASHING_H
#define SPATIAL_HASHING_H

#include "sectors.h"
#include "config.h"
#include "particle.h"

int getCellKey(int cellX, int cellY, SimConfig* config); // Generates a unique cell key based on cell coordinates and configuration
void updateCellKey(Sector* sectors, Particle* particle, SimConfig* config); // Updates the particles cellKey and cellX/Y based on its position and the spatial hashing grid configuration
void updateSpatialLookup(Sector* sectors, Particle particles[], SimConfig* config);

#endif // SPATIAL_HASHING_H