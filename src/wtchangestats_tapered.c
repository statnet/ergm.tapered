#include "ergm_storage.h"
#include "ergm_util.h"
#include "ergm_wtchangestat_operator.h"
#include "ergm_wtchangestats_operator.h"

WtC_CHANGESTAT_FN(c_wttaper_term){
  GET_AUX_STORAGE(StoreWtModelAndStats, storage);
  double *nws = INPUT_PARAM + N_CHANGE_STATS - 1; // Skip the penalty coefficients.

  WtModel *m = storage->m;
  WtChangeStats1(tail, head, weight, nwp, m, edgestate);
  
  memcpy(CHANGE_STAT, m->workspace, (N_CHANGE_STATS-1)*sizeof(double));

  // So only need to change the penalty term statistics
  CHANGE_STAT[N_CHANGE_STATS-1] = 0;
  for(unsigned int k=0; k<m->n_stats; k++){
    double old_diff = storage->stats[k] - nws[k],
           new_diff = old_diff + m->workspace[k];
    CHANGE_STAT[N_CHANGE_STATS-1] += INPUT_PARAM[k]*(new_diff*new_diff - old_diff*old_diff);
  }
}
