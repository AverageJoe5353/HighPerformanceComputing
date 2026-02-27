#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define PI 3.14159265359
#define SIGMA 1.0
#define EPSILON 1.0


int main(int argc, char *argv[]){

    // Parameters
    int N = 100; // Number of Particles
    float ax[N], ay[N], vx[N], vy[N], x[N], y[N]; // Particle Properties
    float vhx[N], vhy[N]; // Half Step Velocities
    float m = 1.0; // Particle Mass

    float L = atof(argv[1]); // Box Length
    float kT = 1.0;

    int i, j; // Loop Variables

    float t=0.0, dt=0.01, tmax=100; // Time Step
    int intsteps = 0; // Integration Steps
    float F_ij=0.0, Fx, Fy, r, dx, dy, theta;
    float F_press=0.0, P=0, P_sum=0, P_sim=0, P_ideal, P_vdw;
    int IS_placed, IS_touching;
    float try_x, try_y;
    float r_cutoff = pow(2.0, 1.0/6.0)*SIGMA;
    
    FILE *outtrj, *outdat;
    outtrj = fopen("gasBox_traj.xyz", "w");
    outdat = fopen("gasBox_pressures.dat", "a");


    // Particle Placement //
    printf("\nInitializing Particles...\n");
    for (i = 0; i < N; i++) {
        IS_placed = 0;
        while (!IS_placed) {
            // Attempt to place particle at random position
            try_x = (float)rand() / (float)RAND_MAX * L - L/2;
            try_y = (float)rand() / (float)RAND_MAX * L - L/2;

            // Test placement attempt for overlap
            IS_touching = 0;
            for (j = 0; j < i; j++) {
                dx = try_x - x[j];
                dy = try_y - y[j];
                r = sqrt(dx * dx + dy * dy);
                
                // If particles are touching, break and try again
                if (r < r_cutoff) {
                    IS_touching = 1;
                    break;
                }
            }
            // If no overlap, place particle
            if (!IS_touching) {
                x[i] = try_x;
                y[i] = try_y;
                IS_placed = 1;
            }
        }

        // Randomize velocity
        theta = (float)rand() / (float)RAND_MAX * 2 * PI;
        vx[i] = sqrt(2 * kT / m) * cos(theta);
        vy[i] = sqrt(2 * kT / m) * sin(theta);

        // Clean Arrays
        ax[i] = 0;
        ay[i] = 0;
        vhx[i] = 0;
        vhy[i] = 0;
    }

    // Make Debug Particles
    // N = 4;
    // x[0] = 2.0;
    // y[0] = 0.0;
    // vx[0] = 0.0;
    // vy[0] = 0.0;

    // x[1] = 0.0;
    // y[1] = 0.0;
    // vx[1] = 0.50;
    // vy[1] = 0.0;

    // x[2] = -2.0;
    // y[2] = 2.0;
    // vx[2] = 0.0;
    // vy[2] = 0.0;

    // x[3] = -2.0;
    // y[3] = 0.0;
    // vx[3] = 0.0;
    // vy[3] = 0.50;



    //// Simulation Loop ////
    printf("Starting Simulation...\n");
    while (t < tmax) {
        // Write positions to trajectory file every 10 steps
        if (intsteps % 10 == 0) {
            fprintf(outtrj, "%i\n", N);
            fprintf(outtrj, "t = %f\n", t);
        }
        
        // Movement Loop
        for (i = 0; i < N; i++) {
            // Write particle positions to traj file every 10 steps
            if (intsteps % 10 == 0) fprintf(outtrj, "a%i %f %f 0\n", i, x[i], y[i]);

            // Half Step Velocity
            vhx[i] = vx[i] + 0.5 * ax[i] * dt;
            vhy[i] = vy[i] + 0.5 * ay[i] * dt;

            // Position Update
            x[i] = x[i] + vhx[i] * dt;
            y[i] = y[i] + vhy[i] * dt;
        }


        // Interaction Loop
        for (i = 0; i < N; i++) {

            // Calculate net force on particle [i]
            Fx = 0.0;
            Fy = 0.0;
            // Get force from all the other particles
            for (j = 0; j < N; j++) {
                if (i != j) {
                    // Calculate distance between particles
                    dx = x[j] - x[i];
                    dy = y[j] - y[i];
                    r = sqrt(dx * dx + dy * dy);

                    if (r < r_cutoff) { // If the particles are within interaction distance (2^(1/6) * sigma)

                        // Calculate force
                        F_ij = -(4.0 * EPSILON) * (1.0/r) * (6.0 * pow(SIGMA/r, 6.0) - 12.0 * pow(SIGMA/r, 12.0));

                        if (F_ij > 1000) {printf("Force over limit at t = %f! F = %f r = %f\n", t, F_ij, r);} // Flag large forces
                        if (F_ij < 0) {printf("Negative Force at t = %f! F = %f r = %f\n", t, F_ij, r);} // Flag negative (attractive) forces

                        Fx -= F_ij * dx / r; // x-component of force
                        Fy -= F_ij * dy / r; // y-component of force
                    }
                }
            }
            // Update Acceleration: a(t+dt)
            ax[i] = Fx / m;
            ay[i] = Fy / m;

            // Update Velocity: v(t+dt)
            vx[i] = vhx[i] + 0.5 * ax[i] * dt;
            vy[i] = vhy[i] + 0.5 * ay[i] * dt;

            // Reflecting Boundary Conditions
            if (x[i] < -L / 2.0 || x[i] > L / 2.0) {
                // Turn around
                vx[i] = -vx[i];
                // Add to pressure force
                F_press += 2*m*fabs(vx[i])/dt;
            }
            if (y[i] < -L / 2.0 || y[i] > L / 2.0) {
                // Turn around
                vy[i] = -vy[i];
                // Add to pressure force
                F_press += 2*m*fabs(vy[i])/dt;
            }
        }


        // Update time
        t += dt;
        intsteps++;

        // Print time every 1000 steps
        if (intsteps % 1000 == 0)
            printf("\r t = %.2f", t);
    }

    // Calculate Pressures
    P_sim = F_press / (intsteps * 4.0 * L);
    P_vdw = N * kT / (L*L - N*2*PI*pow(SIGMA/2,2));
    P_ideal = N * kT / (L*L);

    printf("\nSimulation Complete!\n");
    printf("Num Particles: %i | Box Size: %f\nSimulation Pressure: %f | VDW Pressure: %f | Ideal Gas Pressure: %f\n", N, L, P_sim, P_vdw, P_ideal);
    fprintf(outdat, "%f %f %f %f\n", L, P_sim, P_vdw, P_ideal);
    fclose(outtrj);
    fclose(outdat);
    
    return 0;
}