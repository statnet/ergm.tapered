#include "ergm_storage.h"
#include "ergm_util.h"
#include "ergm_changestat_operator.h"
#include "ergm_changestats_operator.h"

C_CHANGESTAT_FN(c_kurt_term){
  GET_AUX_STORAGE(StoreModelAndStats, storage);
  double *nws = INPUT_PARAM + N_CHANGE_STATS; // Skip the penalty coefficients.

  Model *m = storage->m;
  ChangeStats1(tail, head, nwp, m, edgestate);
  
  memcpy(CHANGE_STAT, m->workspace, (N_CHANGE_STATS)*sizeof(double));

  for(unsigned int k=0; k<m->n_stats; k++){
    double old_diff = storage->stats[k] - nws[k],
           new_diff = old_diff + m->workspace[k];
    double tmp1 = new_diff*new_diff,
           tmp2 = old_diff*old_diff;
    CHANGE_STAT[k           ] = (tmp1 - tmp2);
    CHANGE_STAT[k+m->n_stats] = (tmp1*tmp1 - tmp2*tmp2);
  }
}
