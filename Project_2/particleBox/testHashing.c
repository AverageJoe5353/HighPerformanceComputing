#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

// Demonstation of a faster neighbor Search with Spatial Hashing

#define NUM_PARTICLES 100 // Number of particles
#define L 10            // Size of the 2D box
#define PI 3.14159265358979323846 // Pi
#define R 1.0 // Interaction radius (also cell size for hash grid)
#define KT 1.0 // Boltzmann constant times temperature
#define R_CUTOFF 0.5 // Minimum distance for particle placement to avoid overlaps

// Struct for Particle Information
typedef struct {
    int id;          // Particle ID

    double x, y;     // Position
    double vx, vy;   // Velocity 
    double theta;    // Direction angle

    int cellX;  // Cell coordinates x
    int cellY;  // Cell coordinates y
    int cellKey;     // Cell key for current cell
} Particle;



// Prototypes
void sortParticlesByCellKey(Particle particles[]); // QuickSort Algorithm for particle sort
void updateSpatialLookup(Particle particles[], int startIndices[]); // Updater for neighbor lookup
void updateCellKeys(Particle particles[]); // Update the cell keys for each particle
void placeParticles(Particle particles[]); // Particle placement function to avoid overlaps




int main(){
    // Array of particles
    Particle particles[NUM_PARTICLES]; // Particle List
    int startIndices[NUM_PARTICLES]; // Start Indices for each cell key


    // Seed the random number generator
    srand(time(NULL));

    // Place particles in a grid at 0.5 spacing
    placeParticles(particles);

    

    updateCellKeys(particles); // Update the cell keys for each particle before writing to file

    placeParticles(particles); // Place particles again to randomize positions for the sorted output

    // Write sorted particles to file
    FILE *file1 = fopen("spatialHashGridDemo_unsorted.txt", "w");
    
    if (file1 == NULL) {
        perror("Failed to open file");
        return 1;
    }

    for (int i = 0; i < NUM_PARTICLES; i++) {
        fprintf(file1, "%d %f %f %d %d %d\n", 
            particles[i].id, 
            particles[i].x, particles[i].y, 
            particles[i].cellX, particles[i].cellY, 
            particles[i].cellKey);
    }

    fclose(file1);

    // int n = 0;
    // for (int i = 0; i < L; i++) {
    //     for (int j = 0; j < L; j++) {
    //         particles[n].id = n;
    //         particles[n].x = 0.5 + i;
    //         particles[n].y = 0.5 + j;
    //         n++;
    //     }
    // }

    // Place one extra particle randomly
    // particles[NUM_PARTICLES - 1].id = n + 1;
    // particles[NUM_PARTICLES - 1].x = L * ((double) rand() / RAND_MAX);      
    // particles[NUM_PARTICLES - 1].y = L * ((double) rand() / RAND_MAX);      



    // Update the cell keys for each particle
    updateSpatialLookup(particles, startIndices);
    



    // Write sorted particles to file
    FILE *file = fopen("spatialHashGridDemo_sorted.txt", "w");

    if (file == NULL) {
        perror("Failed to open file");
        return 1;
    }

    for (int i = 0; i < NUM_PARTICLES; i++) {
        fprintf(file, "%d %f %f %d %d %d\n", 
            particles[i].id, 
            particles[i].x, particles[i].y, 
            particles[i].cellX, particles[i].cellY, 
            particles[i].cellKey);
    }

    fclose(file);
    return 0;
}

////////////////////////////////////// FUNCTIONS //////////////////////////////////////

//////////////////////////// SPATIAL HASHING ////////////////////////////
// Discretizing the space into cells and assigning particles to cells.
// A unique key is generated for each cell to store the particles in a hash table.
// Once the particles are sorted by cell key, the start index for each cell key is stored in an array.
// The start index for each group of particles is easily searched in the array via the cell key.

int getCellCoord(float position){
    // Just returns a cell coordinate for a given position
    // Cell size is determined by interaction radius R (Although if I do variable detection sizes, this will need to be changed)
    return (int) (position / R); 
}

int hashCell(int cellX, int cellY){
    // Hash the cell coordinates to get a unique key for the cell
    int primex = 19507;
    int primey = 9737333;
    int numCellY = (int) (L / R); // Number of cells in the y direction (L/R)
    // return (cellX * primex) + (cellY * primey);
    return cellY * numCellY + cellX; // Simple hash function for 2D grid, where numCellY is the number of cells in the y direction (L/R)
}

int getCellKey(int cellHash){
    // Just wraps the cell hash to return a valid index for the cell key
    return cellHash % NUM_PARTICLES;
}

void updateCellKeys(Particle particles[]){
    // Update the cell keys for each particle
    for (int i = 0; i < NUM_PARTICLES; i++) {
        particles[i].cellX = getCellCoord(particles[i].x); // Update the cell coordinates for x
        particles[i].cellY = getCellCoord(particles[i].y); // Update the cell coordinates for y
        // Update the cell key for the particle
        particles[i].cellKey = getCellKey(hashCell(particles[i].cellX, particles[i].cellY));
    }
}

void updateStartIndices(Particle particles[], int startIndices[]){
    // Get Start Indices for cellKeys
    for (int i = 0; i < NUM_PARTICLES; i++) {
        int key = particles[i].cellKey; // Get the cell key of the ith particle

        // If the current key[i] is different from the previous key, i is the start index for that key[i]
        // Special case for the first particle, which is always the start index for the first key.
        if (i == 0 || particles[i - 1].cellKey != key) {
            startIndices[key] = i;
        }
    }
    // startIndices now contains the start index for each cell key group at startIndices[key]
}

// Wrapper function for updating cell keys and start indices
void updateSpatialLookup(Particle particles[], int startIndices[]) {
    updateCellKeys(particles);
    sortParticlesByCellKey(particles);
    updateStartIndices(particles, startIndices);
}
///////////////////////////////////////////////////////////////////////


/////////////////////// QUICKSORT ALGORITHM ///////////////////////
// Just a standard quicksort algorithm to for sorting the particles by cell key.
// By doing this sort the particles are grouped in their array by the cell in which they are located.

void swapParticles(Particle *a, Particle *b) {
    // Swaps the array placement for two particles
    Particle temp = *a; // Store the value of a
    *a = *b; // Move the value of b to a's spot
    *b = temp; // Move the value of the original a to b's spot
}

int partition(Particle particles[], int low, int high) {
    // Sorting particles by cell keys
    // Split list around a pivot point
    // Move all smaller values to the left of the pivot, larger to the right

    // Choose a pivot point
    int pivot = particles[high].cellKey; // Pivot around the cell key of the last particle

    // Initialize the index of the smaller element, i
    int i = low - 1;

    // Loop through the list of particles, comparing cell keys to the pivot
    for (int j = low; j < high; j++) {
        if (particles[j].cellKey < pivot) { // If the cell key of particle[j] is less than the pivot value
            i++; // Increment i (so it places above low)
            swapParticles(&particles[i], &particles[j]); // swap the current particle with the particle at index i
        }
    }

    // After the loop, all particles with cell keys less than the pivot are to the left of i
    // So swap the particle at i + 1 with the particle at the pivot
    swapParticles(&particles[i + 1], &particles[high]);

    // Return the index of the new pivot
    return (i + 1);
}

void quickSort(Particle particles[], int low, int high) {
    // Recursive function to sort particles by cell key
    if (low < high) {
        // Partition the array, pi is the returned pivot index from partition
        int pi = partition(particles, low, high);

        // Recursively sort the subarrays
        quickSort(particles, low, pi - 1);
        quickSort(particles, pi + 1, high);

        // Continues recursively until the low and pi-1 are equal (or pi+1 and high)
        // At which point the array is sorted
    }
}

void sortParticlesByCellKey(Particle particles[]) {
    // Wrapper function to sort particles by cell key
    quickSort(particles, 0, NUM_PARTICLES - 1);
}

///////////////////////////////////////////////////////////////////////

void placeParticles(Particle particles[]) {
    double halfBox = L / 2.0;
    for (int i = 0; i < NUM_PARTICLES; i++) {
        int overlap;
        do {
            overlap = 0;
            // Randomly place particle within the box
            particles[i].x = ((double)rand() / RAND_MAX) * L;
            particles[i].y = ((double)rand() / RAND_MAX) * L;

            // Check for overlaps with previously placed particles
            for (int j = 0; j < i; j++) {
                double dx = particles[i].x - particles[j].x;
                double dy = particles[i].y - particles[j].y;
                double r2 = dx * dx + dy * dy;
                if (r2 < R_CUTOFF * R_CUTOFF) {
                    overlap = 1; // Overlap detected, need to reposition
                    break;
                }
            }
        } while (overlap); // Repeat until no overlaps

        double ran1=(float)rand()/RAND_MAX;
        double ran2=(float)rand()/RAND_MAX;

        particles[i].vx = sqrt(KT/(1*NUM_PARTICLES))*sqrt(-2*log(ran1))*cos(2*PI*ran2);
        particles[i].vy = sqrt(KT/(1*NUM_PARTICLES))*sqrt(-2*log(ran1))*sin(2*PI*ran2);

    }
}