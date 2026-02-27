#include "quicksort.h"

#include "particle.h"

/////////////////////// QUICKSORT ALGORITHM ///////////////////////
// Just a standard quicksort algorithm to for sorting the particles and obstacles by cell key.
// By doing this sort the particles are grouped in their array by the cell in which they are located.
// Modified from: https://www.geeksforgeeks.org/quick-sort-algorithm/



void swap(Particle *a, Particle *b) {
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
            swap(&particles[i], &particles[j]); // swap the current particle with the particle at index i
        }
    }

    // After the loop, all particles with cell keys less than the pivot are to the left of i
    // So swap the particle at i + 1 with the particle at the pivot
    swap(&particles[i + 1], &particles[high]);

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


///////////////////////////////////////////////////////////////////////

