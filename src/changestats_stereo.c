#include "ergm_storage.h"
#include "ergm_util.h"
#include "ergm_changestat_operator.h"
#include "ergm_changestats_operator.h"

C_CHANGESTAT_FN(c_stereo_term){
  GET_AUX_STORAGE(StoreModelAndStats, storage);
  double *nws = INPUT_PARAM + N_CHANGE_STATS - 1; // Skip the penalty coefficients.
  double new_s, cur_s;

  Model *m = storage->m;
  ChangeStats1(tail, head, nwp, m, edgeflag);
  
  memcpy(CHANGE_STAT, m->workspace, (N_CHANGE_STATS-1)*sizeof(double));

  new_s = INPUT_PARAM[0]*INPUT_PARAM[0];
  cur_s = INPUT_PARAM[0]*INPUT_PARAM[0];
  for(unsigned int k=0; k<m->n_stats; k++){
    double old_diff = storage->stats[k] - nws[k],
           new_diff = old_diff + m->workspace[k];
    new_s += new_diff*new_diff;
    cur_s += old_diff*old_diff;
  }

  CHANGE_STAT[N_CHANGE_STATS-1] = 2.0*log(cur_s) - 2.0*log(new_s);
}
