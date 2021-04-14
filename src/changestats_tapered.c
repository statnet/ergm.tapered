#include "ergm_storage.h"
#include "ergm_util.h"
#include "ergm_changestat_operator.h"
#include "ergm_changestats_operator.h"

C_CHANGESTAT_FN(c_taper_term){
  GET_AUX_STORAGE(StoreModelAndStats, storage);
  double *nws = INPUT_PARAM + N_CHANGE_STATS - 1; // Skip the beta penalty coefficients.
     // The penalty term adds an extra change statistic (for the -1)

//Rprintf("INPUT_PARAM %d N_CHANGE_STATS %d\n", INPUT_PARAM, N_CHANGE_STATS);
//Rprintf(" %f %f %f  %f %f %f",nws[0],nws[1],nws[2],nws[3],nws[4],nws[5]);

  Model *m = storage->m;
  ChangeStats1(tail, head, nwp, m, edgeflag);
  
  // Copy the change statistics to CHANGE_STAT
  memcpy(CHANGE_STAT, m->workspace, (N_CHANGE_STATS-1)*sizeof(double));

  // So only need to change the penalty term statistics
  CHANGE_STAT[N_CHANGE_STATS-1] = 0.0;
  for(unsigned int k=0; k<m->n_stats; k++){
    double old_diff = storage->stats[k] - nws[k],
           new_diff = old_diff + m->workspace[k];
//  Rprintf("k %d nws[k] %f old_diff %f new_diff %f",k,nws[k],old_diff,new_diff);
//if(k==6) Rprintf("k %d i[k] %f orig_nws[k] %f cur_nws[k] %f new_nws[k] %f old_diff %f new_diff %f\n",k,INPUT_PARAM[k],nws[k],storage->stats[k],storage->stats[k]+ m->workspace[k],old_diff,new_diff);
//  if(k==3 && (storage->stats[k] + m->workspace[k] < 0.0)) Rprintf("k %d i[k] %f orig_nws[k] %f cur_nws[k] %f new_nws[k] %f old_diff %f new_diff %f\n",k,INPUT_PARAM[k],nws[k],storage->stats[k],storage->stats[k]+ m->workspace[k],old_diff,new_diff);
    CHANGE_STAT[N_CHANGE_STATS-1] += INPUT_PARAM[k]*(new_diff*new_diff - old_diff*old_diff);
  }
//  Rprintf("\n");
}
