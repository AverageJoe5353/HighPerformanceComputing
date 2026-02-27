#ifndef SECTORS_H
#define SECTORS_H

typedef struct {
    int key;
    int startIndex;
    int population;
    int count;
} Sector;

Sector makeSectors();

#endif // SECTORS_H