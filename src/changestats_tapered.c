#include "ergm_storage.h"
#include "ergm_util.h"
#include "ergm_changestat_operator.h"

typedef struct{Model *m; double *stats;} StoreModelAndStats;

I_CHANGESTAT_FN(i_taper_term){
  double *inputs = INPUT_PARAM + N_CHANGE_STATS - 1; // Skip the penalty coefficients.
  ALLOC_STORAGE(1, StoreModelAndStats, storage);

  // Unpack the submodel.
  Model *m = storage->m = ModelInitialize(getListElement(mtp->R, "submodel"), NULL,  nwp, FALSE);

  storage->stats = Calloc(m->n_stats, double);
  memcpy(storage->stats, inputs, m->n_stats*sizeof(double));

  SummStats(0, NULL, NULL, nwp, m);
  addonto(storage->stats, m->workspace, m->n_stats);

  DELETE_IF_UNUSED_IN_SUBMODEL(z_func, m);
}

C_CHANGESTAT_FN(c_taper_term){
  GET_STORAGE(StoreModelAndStats, storage);
  Model *m = storage->m;

  // TODO: cache the result for use in u_taper_term.
  ChangeStats1(tail, head, nwp, m, edgeflag);
  
  memcpy(CHANGE_STAT, m->workspace, (N_CHANGE_STATS-1)*sizeof(double));

  CHANGE_STAT[N_CHANGE_STATS-1] = 0;
  for(unsigned int k=0; k<m->n_stats; k++){
    double new_diff = storage->stats[k]+m->workspace[k];
    CHANGE_STAT[N_CHANGE_STATS-1] += INPUT_PARAM[k]*(new_diff*new_diff - storage->stats[k]*storage->stats[k]);
  }
}

U_CHANGESTAT_FN(u_taper_term){
  GET_STORAGE(StoreModelAndStats, storage);
  Model *m = storage->m;

  // TODO: Use the result cached by the c_ function.
  ChangeStats1(tail, head, nwp, m, edgeflag);
  addonto(storage->stats, m->workspace, m->n_stats);
}

F_CHANGESTAT_FN(f_taper_term){
  GET_STORAGE(StoreModelAndStats, storage);

  Free(storage->stats);
  ModelDestroy(nwp, storage->m);
}
