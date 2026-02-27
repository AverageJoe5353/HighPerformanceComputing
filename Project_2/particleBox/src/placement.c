#include "placement.h"

#include "config.h"
#include "particle.h"

#include <stdlib.h>
#include <math.h>

double dx, dy, r;

int i, j;

#define PI 3.14159265359f

void placeParticles(Particle particles[], SimConfig* config) {
    double halfBox = config->boxLength / 2.0;
    for (int i = 0; i < config->numParticles; i++) {
        int overlap;
        do {
            overlap = 0;
            // Randomly place particle within the box
            particles[i].x = ((double)rand() / RAND_MAX) * config->boxLength - halfBox;
            particles[i].y = ((double)rand() / RAND_MAX) * config->boxLength - halfBox;

            // Check for overlaps with previously placed particles
            for (int j = 0; j < i; j++) {
                double dx = particles[i].x - particles[j].x;
                double dy = particles[i].y - particles[j].y;
                double r2 = dx * dx + dy * dy;
                if (r2 < config->r_cutoff * config->r_cutoff) {
                    overlap = 1; // Overlap detected, need to reposition
                    break;
                }
            }
        } while (overlap); // Repeat until no overlaps

        double ran1=(float)rand()/RAND_MAX;
        double ran2=(float)rand()/RAND_MAX;

        particles[i].vx = sqrt(config->kT/(particles[i].mass*config->numParticles))*sqrt(-2*log(ran1))*cos(2*PI*ran2);
        particles[i].vy = sqrt(config->kT/(particles[i].mass*config->numParticles))*sqrt(-2*log(ran1))*sin(2*PI*ran2);

        particles[i].ax = 0.0;
        particles[i].ay = 0.0;
    }
}