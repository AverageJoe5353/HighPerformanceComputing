/*
This is a randomized placement algorithm for particles that avoids overlaps.

It is designed to be a callable function in main(), but you can also drop it directly into main().

If you are adding this as a callable function, copy the whole thing and paste it into the .c file where you want it.
If you are adding this directly into main(), copy the contents of the function without the `void placeParticles(...) {` and `}` lines and paste it where you want the placement to happen.
*/


// Required Libraries (make sure these are included at the top of main.c or whatever file you are putting this in)
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Required Constants (make sure these are defined at the top of main.c or somewhere this function can see them) 
#define PI 3.14159265359f
#define SIGMA 1.0
#define N 1000  // Number of Particles
#define L 100.0 // Box Length
#define kT 1.0 // Temperature

// Particle Placement Algorithm
void placeParticles(float *x, float *y, float *vx, float *vy, float *ax, float *ay, float *vhx, float *vhy, float *m) {
    /* Pass the arrays for each value as arguments */
    /* The arrays should be initialized before calling this function like this: double x[N], y[N], and so on... */
    /* Parameters like N, L, SIGMA, and kT are defined as macros above, but you might have to tweak some stuff to get the function to work */


    int i, j;
    float try_x, try_y, dx, dy, r, theta;
    int IS_placed, IS_touching;
    float r_cutoff = pow(2, 1.0/6.0) * SIGMA; // Overlap cutoff distance

    // Loop through each particle and attempt to place it without overlap
    for (i = 0; i < N; i++) {
        // Set initial values for this particle
        ax[i] = 0.0;
        ay[i] = 0.0;
        vhx[i] = 0.0;
        vhy[i] = 0.0;
        m[i] = 1.0; // Set unit mass for all particles - comment this out if you have masses set already


        IS_placed = 0;
        while (!IS_placed) {
            // Attempt to place particle at random position
            try_x = (float)rand() / (float)RAND_MAX * L - L/2;
            try_y = (float)rand() / (float)RAND_MAX * L - L/2;

            // Test placement attempt for overlap
            IS_touching = 0;
            for (j = 0; j < i; j++) { // Loop through already placed particles
                dx = try_x - x[j];
                dy = try_y - y[j];
                r = sqrt(dx * dx + dy * dy);
                
                // If particles are touching, break and try again
                if (r < r_cutoff) {
                    IS_touching = 1;
                    break; // Breaks the j loop at first detection of overlap
                }
            }
            // If no overlap, place particle
            if (!IS_touching) {
                x[i] = try_x;
                y[i] = try_y;
                IS_placed = 1; // Yay!
            }
        }

        // Randomize velocity
        theta = (float)rand() / (float)RAND_MAX * 2 * PI;
        vx[i] = sqrt(2 * kT / m[i]) * cos(theta);
        vy[i] = sqrt(2 * kT / m[i]) * sin(theta);

        
    }
}