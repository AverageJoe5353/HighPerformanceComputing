#include "forces.h"

#include "config.h"
#include "particle.h"

#include <math.h>
#include <stdio.h>

static double dx, dy, dz, r2, r_cutoff2, r6, r12, sigma6, sigma12, F_ij;
#define PI 3.14159265358979323846

void applyForce_LJ(Particle* p1, Particle* p2, SimConfig* config) {
    dx = p1->x - p2->x;
    dy = p1->y - p2->y;
    dz = p1->z - p2->z;

    r2 = dx * dx + dy * dy + dz * dz;
    r_cutoff2 = config->r_cutoff * config->r_cutoff;

    if (r2 < r_cutoff2) {
        // printf("Collision detected between Particle %d and Particle %d\n", p1->ID, p2->ID);
        r6 = r2 * r2 * r2;
        r12 = r6 * r6;
        sigma6 = pow(config->sigma, 6);
        sigma12 = pow(config->sigma, 12);
        F_ij = -(4.0 * config->epsilon) * (1.0/sqrt(r2)) * ((6.0 * sigma6 / r6) - (12.0 * sigma12 / r12));

        // Update net forces on each particle
        p1->Fx += F_ij * dx / sqrt(r2);
        p1->Fy += F_ij * dy / sqrt(r2);
        p1->Fz += F_ij * dz / sqrt(r2);
    } else {
        // No force if particles are beyond cutoff distance
        F_ij = 0.0;
    }
}

void applyReflectingBoundary(Particle* p, SimConfig* config) {
    double halfBox = config->boxLength / 2.0;

    while (p->x < -halfBox || p->x > halfBox) {
        // printf("Particle %d escaped boundary in x direction at position %f\n", p->ID, p->x);
        if (p->x < -halfBox) {
            p->x = -2*halfBox - (p->x); // Reflect position
            p->vx = -p->vx; // Reverse velocity
        } else if (p->x > halfBox) {
            p->x =  2*halfBox - (p->x); // Reflect position
            p->vx = -p->vx; // Reverse velocity
        }
    }

    while (p->y < -halfBox || p->y > halfBox) {
        // printf("Particle %d escaped boundary in y direction at position %f\n", p->ID, p->y);
        if (p->y < -halfBox) {
            p->y = -2*halfBox - (p->y); // Reflect position
            p->vy = -p->vy; // Reverse velocity
        } else if (p->y > halfBox) {
            p->y =  2*halfBox - (p->y); // Reflect position
            p->vy = -p->vy; // Reverse velocity
        }
    }

    while (p->z < -halfBox || p->z > halfBox) {
        // printf("Particle %d escaped boundary in z direction at position %f\n", p->ID, p->z);
        if (p->z < -halfBox) {
            p->z = -2*halfBox - (p->z); // Reflect position
            p->vz = -p->vz; // Reverse velocity
        } else if (p->z > halfBox) {
            p->z =  2*halfBox - (p->z); // Reflect position
            p->vz = -p->vz; // Reverse velocity
        }
    }
}

double d1(float x, float y, float z, float a){
    // Distance to left edge of pore located at (-a,0)
    double rho = sqrt(x*x + y*y);
    return sqrt((rho + a)*(rho + a) + z*z);
}
double d2(float x, float y, float z, float a){
    // Distance to right edge of pore located at (a,0)
    double rho = sqrt(x*x + y*y);
    return sqrt((rho - a)*(rho - a) + z*z);
}

void applyForce_nanoporeEfield(Particle* p, SimConfig* config) {
    double dist1 = d1(p->x, p->y, p->z, config->pore_a);
    double dist2 = d2(p->x, p->y, p->z, config->pore_a);
    double r = sqrtl(0.5*(dist1+dist2)-(config->pore_a*config->pore_a));
    double Exy = -(config->V/PI) * ((dist1-dist2)/(dist1+dist2)) * (1/(dist1*dist2)) * sqrtl(((dist1+dist2)/2)*((dist1+dist2)/2)-(config->pore_a*config->pore_a));
    double Ez = -(config->V/PI) * (1/(dist1*dist2)) * sqrtl((config->pore_a*config->pore_a)-((dist1-dist2)/2)*((dist1-dist2)/2));
    
    p->Fx += p->charge * Exy * (p->x / r);
    p->Fy += p->charge * Exy * (p->y / r);
    p->Fz += p->charge * Ez;
}