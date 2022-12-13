#ifndef __SE_FKG_H
#define __SE_FKG_H

#include "SE_fgg.h"

#ifdef THREE_PERIODIC
#define __FKG_EXPA kaiser_expansion_3p
#ifdef FORCE
#define __FKG_EXPA_FORCE kaiser_expansion_3p_force
#endif
#define __FKG_INDEX kaiser_index_3p
#endif

#ifdef TWO_PERIODIC
#define __FKG_EXPA kaiser_expansion_2p
#ifdef FORCE
#define __FKG_EXPA_FORCE kaiser_expansion_2p_force
#endif
#define __FKG_INDEX kaiser_index_2p
#endif

#ifdef ONE_PERIODIC
#define __FKG_EXPA kaiser_expansion_1p
#ifdef FORCE
#define __FKG_EXPA_FORCE kaiser_expansion_1p_force
#endif
#define __FKG_INDEX kaiser_index_1p
#endif

// kaiser signatures
void SE_FKG_allocate_workspace(SE_FGG_work*, const SE_FGG_params*, int);
void SE_FKG_expand_all(SE_FGG_work*, const SE_state* , const SE_FGG_params*);
void SE_FKG_expand_all_force(SE_FGG_work*, const SE_state* , const SE_FGG_params*);
void SE_FKG_grid(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FKG_int(double*, const SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FKG_int_force(double*, const SE_FGG_work*, const SE_state*, const SE_FGG_params*);

void SE_FKG_int_split_SSE_dispatch(double*, const SE_FGG_work*, const SE_FGG_params*);
void SE_FKG_grid_split_SSE_dispatch(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FKG_int_split_AVX_dispatch(double*, const SE_FGG_work*, const SE_FGG_params*);
void SE_FKG_grid_split_AVX_dispatch(SE_FGG_work*, const SE_state*, const SE_FGG_params*);
void SE_FKG_int_split_force(double*, SE_state*, const SE_FGG_work*, const SE_FGG_params*);

#endif
