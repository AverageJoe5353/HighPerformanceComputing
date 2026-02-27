#include "sectors.h"

Sector makeSectors() {
    Sector s;
    s.key = -1;
    s.startIndex = -1;
    s.population = 0;
    s.count = 0;
    return s;
}