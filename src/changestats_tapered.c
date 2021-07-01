#include "ergm_storage.h"
#include "ergm_util.h"
#include "ergm_changestat_operator.h"
#include "ergm_changestats_operator.h"

C_CHANGESTAT_FN(c_taper_term){
  GET_AUX_STORAGE(StoreModelAndStats, storage);
  double *nws = INPUT_PARAM + N_CHANGE_STATS - 1; // Skip the beta penalty coefficients.
     // The penalty term adds an extra change statistic (for the -1)

  Model *m = storage->m;
  ChangeStats1(tail, head, nwp, m, edgestate);
  
  // Copy the change statistics to CHANGE_STAT
  memcpy(CHANGE_STAT, m->workspace, (N_CHANGE_STATS-1)*sizeof(double));

  // So only need to change the penalty term statistics
  CHANGE_STAT[N_CHANGE_STATS-1] = 0.0;
  for(unsigned int k=0; k<m->n_stats; k++){
    double old_diff = storage->stats[k] - nws[k],
           new_diff = old_diff + m->workspace[k];
    CHANGE_STAT[N_CHANGE_STATS-1] += INPUT_PARAM[k]*(new_diff*new_diff - old_diff*old_diff);
  }
}
