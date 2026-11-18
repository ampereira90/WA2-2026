/*
 * ACA - Work Assignment 2 - variant "box27": stencil coefficients
 */
#include "box27.h"

#include <cstdlib>

/* weight of neighbour (di,dj,dk) in the 27-point operator (see box27.h) */
double stencil_weight(int di, int dj, int dk)
{
    int dist = std::abs(di) + std::abs(dj) + std::abs(dk);
    switch (dist) {
        case 0:  return 128.0 / 30.0;   /* centre */
        case 1:  return  14.0 / 30.0;   /* face   */
        case 2:  return   3.0 / 30.0;   /* edge   */
        case 3:  return   1.0 / 30.0;   /* corner */
        default: return 0.0;
    }
}
