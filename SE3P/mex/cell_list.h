#ifndef __CELL_LIST_H__
#define __CELL_LIST_H__

#include <stdio.h>

void build_cell_list(
                     // Input
                     const double* restrict x,
                     const int N,
                     const double* restrict box,
                     const double rc,
                     // Output
                     double rn[3],
                     int ncell[3],
                     int* restrict *cell_list_p,
                     int* restrict *cell_idx_p
                     );

#endif
