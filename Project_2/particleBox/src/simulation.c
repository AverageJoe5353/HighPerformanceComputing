#include "simulation.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "spatialHashing.h"
#include "forces.h"
#include "sectors.h"
#include "config.h"
#include "particle.h"


static int p,j;
int t;
static double drag, noise;


void runSimulation(Particle particles[], SimConfig* config) {
    // Create sector list
    Sector sectors[config->numCells];
    for (int i = 0; i < config->numCells; i++) {
        sectors[i] = makeSectors();
    }

    drag = config->gamma;
    noise = sqrt(2 * config->kT * config->gamma / config->timeStep);

    FILE *outt,*outl,*outstats,*outpressure;
	outt=fopen("langevin2D_traj.xyz","w");
    outl=fopen("langevin2D_traj.lammpstrj","w");
    outpressure=fopen("langevin2D_pressure.dat","w");
	outstats=fopen("langevin2D_stats.txt","w");
    
    while (t < config->maxStep) {
        // if(t%config->recordInterval==0 && config->recordToggle){
        //     // Write to traj header
        //     // printf("\rStep %3d / %4d", t, (int)config->maxStep);
		// 	fprintf(outt,"%i\n",config->numParticles);
		// 	fprintf(outt,"%f\n",t * config->timeStep);
		// }

        if (t%config->recordInterval==0 && config->recordToggle){
            // Write to lammpstrj header
            fprintf(outl,"ITEM: TIMESTEP\n");
            fprintf(outl,"%i\n",t);
            fprintf(outl,"ITEM: NUMBER OF ATOMS\n");
            fprintf(outl,"%i\n",config->numParticles);
            fprintf(outl,"ITEM: BOX BOUNDS pp pp pp\n");
            fprintf(outl,"0 %f\n",config->boxLength);
            fprintf(outl,"0 %f\n",config->boxLength);
            fprintf(outl,"0 %f\n",config->boxLength);
            fprintf(outl,"ITEM: ATOMS id type fragment x y z vx vy vz\n");
        }

        debugrep("Starting movement phase...");
        // Reset Sectors
        for (int i = 0; i < config->numCells; i++) {
            sectors[i].population = 0;
            sectors[i].count = 0;
        }
        // Action Phase (Position Update)
        for (p = 0; p < config->numParticles; p++) {

            // if(t%config->recordInterval==0 && config->recordToggle){
			// 	// Write to traj
            //     fprintf(outt,"B%i %f %f %f\n",particles[p].ID,particles[p].x,particles[p].y,particles[p].z);
            // }
            
            if (t%config->recordInterval==0 && config->recordToggle){
                // Write to lammpstrj
                fprintf(outl,"%i A %i %f %f %f %i %i %i\n",particles[p].ID, particles[p].cellKey, particles[p].x, particles[p].y, particles[p].z, particles[p].cellX, particles[p].cellY*10, particles[p].cellKey);
            }

            // Update Half-Step Velocity
            debugstage("Updating half-step velocities...");
            particles[p].vhx =  particles[p].vx + 0.5 * particles[p].ax * config->timeStep;
            particles[p].vhy =  particles[p].vy + 0.5 * particles[p].ay * config->timeStep;

            // Update Positions
            debugstage("Updating positions...");
            particles[p].x += particles[p].vhx * config->timeStep;
            particles[p].y += particles[p].vhy * config->timeStep;

            // Apply Boundary Conditions
            debugstage("Applying boundary conditions...");
            applyReflectingBoundary(&particles[p], config);

            

            // Update Cell Key
            debugstage("Updating cell keys...");    
            updateCellKey(sectors, &particles[p], config);
        }


        debugrep("Updating spatial lookup...");
        updateSpatialLookup(sectors, particles, config); // O(n log n) spatial hashing update (worst case O(N^2))

        // Interaction Phase
        debugrep("Starting interaction phase...");
        for (p = 0; p < config->numParticles; p++) {
            particles[p].Fx = 0.0;
            particles[p].Fy = 0.0;
            // Collisions with Neighbouring Particles (LJ Force)
            for (j = 0; j < 9; j++) { // Loop over neighbouring cells
                // printf("Checking neighbouring cell %d for Particle %d\n", particles[p].cellSearchKey[j], particles[p].ID);
                int searchKey = particles[p].cellSearchKey[j];
                
                if (searchKey != -1) { // If the neighbouring cell exists AND has particles in it
                    int startIndex = sectors[searchKey].startIndex;
                    if (startIndex != -1) { // If there are particles in the neighbouring cell
                        for (int k = startIndex; k < config->numParticles && particles[k].cellKey == searchKey; k++) { // Loop over particles in neighbouring cell
                            if (particles[p].ID != particles[k].ID){
                                applyForce_LJ(&particles[p], &particles[k], config);
                                // printf("Checking Particle %d against Particle %d\n", particles[p].ID, particles[k].ID);
                            }
                        }
                    }
                }
            }

            // External Forces (e.g. Gravity, Electric Field) if wanted...
            // applyForce_nanoporeEfield(&particles[p], config);
        }




        // Reaction Phase (Acceleration and Velocity Update)
        debugrep("Starting reaction phase...");
        for (p = 0; p < config->numParticles; p++) {
            
            // Update Accelerations
            particles[p].ax = -drag*particles[p].vhx + noise*((double)rand() / RAND_MAX - 0.5) + particles[p].Fx;
            particles[p].ay = -drag*particles[p].vhy + noise*((double)rand() / RAND_MAX - 0.5) + particles[p].Fy;

            // Update Full-Step Velocity
            particles[p].vx =  particles[p].vx + 0.5 * particles[p].ax * config->timeStep;
            particles[p].vy =  particles[p].vy + 0.5 * particles[p].ay * config->timeStep;

        }

    t++;
    printf("\rStep %3d / %4d", t, (int)config->maxStep);
    }
    printf("\n");
    printf("Total steps: %d\n", t);
}