#include "ergm_storage.h"
#include "ergm_util.h"
#include "ergm_changestat_operator.h"
#include "ergm_changestats_operator.h"

C_CHANGESTAT_FN(c_var_term){
  GET_AUX_STORAGE(StoreModelAndStats, storage);
  double *nws = INPUT_PARAM; // The network statistics are up first

  Model *m = storage->m;
  ChangeStats1(tail, head, nwp, m, edgeflag);
  
  // Copy the change statistics to CHANGE_STAT
  memcpy(CHANGE_STAT, m->workspace, (N_CHANGE_STATS)*sizeof(double));

  for(unsigned int k=0; k<m->n_stats; k++){
    double old_diff = storage->stats[k] - nws[k],
           new_diff = old_diff + m->workspace[k];
    CHANGE_STAT[k] = INPUT_PARAM[k]*(new_diff*new_diff - old_diff*old_diff);
  }
}
