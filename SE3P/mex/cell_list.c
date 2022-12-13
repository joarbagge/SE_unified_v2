#include "cell_list.h"
#include <math.h>
#include "mex_compat.h"
#define MIN(a, b) ((a > b) ? b : a)

/* =============================================================================
 * ==== BUILD CELL LIST ========================================================
 *
 * Divides the computational box into cells, and creates a list
 * of the particles in each cell. All cells have the same size,
 * and all their side lengths are guaranteed to be at least rc.
 * It is assumed that rc <= min(box). Otherwise the function will
 * abort, and the targets of the cell_list_p and cell_idx_p pointers
 * will both be set to NULL.
 *
 * The cells will not necessarily be cubic (although they are if
 * the box itself is a cube).
 *
 * Input:
 * - x: positions of the particles in the box
 * - N: number of particles
 * - box: size of computational box [L1,L2,L3]
 * - rc: the cut-off radius
 *
 * Output:
 * - rn: array of side lengths of each cell [rn1,rn2,rn3]
 * - ncell: array of number of cells in each direction [nc1,nc2,nc3]
 * - cell_list_p: pointer to an array (of size N) containing a
 *   permutation of particle indices
 * - cell_idx_p: pointer to an array (of size prod(ncell)+1)
 *   containing the cumulative count of particles in each cell
 *
 * In the code, `i` is the index of a cell (sorted by their
 * location in the box), and `n` is the index of a particle
 * (by their order in the `x` array). Furthermore, let `pidx`
 * be the index of a particle, sorted by the cell the particle
 * belongs to. Then,
 * - Particles in cell i has indices pidx such that cell_idx[i] <= pidx < cell_idx[i]
 * - To convert to index n, use n = cell_list[pidx]
 */
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
                     )
{
    // Output (will return pointers to these pointers)
    int* restrict cell_list;
    int* restrict cell_idx;
    // Intermediates (could do this with fewer variables, but this is clear)
    int* restrict point_cell_map; // which cell each particle is in
    int* restrict cell_count; // number of particles in each cell

    // Setup cell partitioning
    double boxmin = MIN(box[0], box[1]);
    boxmin = MIN(boxmin, box[2]);
    if (rc > boxmin)
    {
        *cell_list_p = NULL;
        *cell_idx_p = NULL;
        return;
    }
    for (int j=0; j<3; j++)
    {
        ncell[j] = floor(box[j] / rc);
        rn[j] = box[j] / ncell[j];
    }
    int ncell_tot = ncell[0]*ncell[1]*ncell[2];

    // Allocate arrays
    cell_list = __MALLOC(N*sizeof(int));
    cell_idx = __CALLOC(ncell_tot+1, sizeof(int));
    point_cell_map = __MALLOC(N*sizeof(int));
    cell_count = __CALLOC(ncell_tot, sizeof(int));

    // Build list in two sweeps.
    // First sweep, fill point_cell_map and cell_idx
    int icell[3]; // index of a cell in each direction
    for (int n=0; n<N; n++)
    {
        // Find icell of the cell that particle n belongs to
        for (int j=0; j<3; j++)
        {
            icell[j] = x[n*3+j]/rn[j];
        }
        // Find index i of this cell
        int i =
            icell[0] +
            icell[1]*ncell[0] +
            icell[2]*ncell[1]*ncell[0];
        point_cell_map[n] = i; // store index of this cell
        cell_idx[i+1]++; // count number of particles in this cell
    }
    // Make the count cumulative
    cell_idx[0] = 0;
    for (int i=0; i<ncell_tot; i++)
    {
        cell_idx[i+1] += cell_idx[i];
    }

    // Second sweep, fill cell_list (and cell_count)
    for (int n=0; n<N; n++)
    {
        int i = point_cell_map[n]; // index of the cell this particle belongs to
        int pidx = cell_idx[i] + cell_count[i];
        cell_list[pidx] = n;
        cell_count[i]++;
    }

    // Free intermediate variables
    __FREE(point_cell_map);
    __FREE(cell_count);

    // Output
    *cell_list_p = cell_list;
    *cell_idx_p = cell_idx;
}
