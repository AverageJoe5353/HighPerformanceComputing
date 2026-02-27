#ifndef QUICKSORT_H
#define QUICKSORT_H

#include "config.h"
#include "particle.h"

void swap(Particle *a, Particle *b);
int partition(Particle particles[], int low, int high);
void quickSort(Particle particles[], int low, int high);

#endif // QUICKSORT_H