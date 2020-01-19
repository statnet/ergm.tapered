#include "ergm_changestat_operator.h"

typedef struct{Model *m; double *stats;} StoreModelAndStats;

I_CHANGESTAT_FN(i_m4_term){
  double *inputs = INPUT_PARAM + N_CHANGE_STATS; // Skip the penalty coefficients.
  ALLOC_STORAGE(1, StoreModelAndStats, storage);

  // Unpack the submodel.
  Model *m = storage->m = unpack_Model_as_double(&inputs);
  // Initialize empty network.
  Network *nw = NetworkInitialize(NULL, NULL, 0, N_NODES, DIRECTED, BIPARTITE, 0, 0, NULL);

  // Initialize storage for submodel terms 
  InitStats(nw, m);

  storage->stats = Calloc(m->n_stats, double);
  memcpy(storage->stats, inputs, m->n_stats*sizeof(double));

  // Evaluate the initial difference between the current nwp and the m vector (the slow way).
  EXEC_THROUGH_NET_EDGES(tail, head, e, {
      ChangeStats(1, &tail, &head, nw, m);
      for(unsigned int k=0; k<m->n_stats; k++)
	storage->stats[k] += m->workspace[k];
      Rboolean edgeflag = IS_OUTEDGE((tail), (head), (nw));
      UPDATE_STORAGE(tail, head, nw, m, NULL, edgeflag);
      ToggleEdge(tail, head, nw);
    });
  
  // Note that nw is now identitical to nwp.
  NetworkDestroy(nw);
}

D_CHANGESTAT_FN(d_m4_term){
  GET_STORAGE(StoreModelAndStats, storage);
  Model *m = storage->m;

  // TODO: cache the result for use in u_m4_term.
  ChangeStats(ntoggles, tails, heads, nwp, m);
  
  memcpy(CHANGE_STAT, m->workspace, (N_CHANGE_STATS)*sizeof(double));

  int nstat = m->n_stats;

//Rprintf("nstat %d\n",nstat);
  for(unsigned int k=0; k<nstat; k++){
    double new_diff = storage->stats[k]+m->workspace[k];
    double tmp1 = new_diff*new_diff;
    double tmp2 = storage->stats[k]*storage->stats[k];
    // Note that this is *negative* the fourth moment
    CHANGE_STAT[k] = INPUT_PARAM[k]*(tmp2*tmp2 - tmp1*tmp1);
//Rprintf("k %d nstat %d tmp1 %f tmp2 %f\n",k,nstat,tmp1,tmp2);
//Rprintf("k %d nstat %d C1 %f C2 %f\n",k,nstat,(tmp1 - tmp2),(tmp1*tmp1 - tmp2*tmp2));
//Rprintf("k %d nstat %d C1 %f C2 %f\n",k,nstat,(tmp1 - tmp2),(tmp1*tmp1 - tmp2*tmp2));
//Rprintf("k %d nstat %d C1 %f C2 %f\n",k,nstat,INPUT_PARAM[k      ]*(tmp1 - tmp2),INPUT_PARAM[k+nstat]*(tmp1*tmp1 - tmp2*tmp2));
  }
}

U_CHANGESTAT_FN(u_m4_term){
  GET_STORAGE(StoreModelAndStats, storage);
  Model *m = storage->m;

  // TODO: Use the result cached by the d_ function.
  ChangeStats(1, &tail, &head, nwp, m);
  for(unsigned int k=0; k<m->n_stats; k++)
    storage->stats[k] += m->workspace[k];

  UPDATE_STORAGE(tail, head, nwp, m, NULL, edgeflag);
}

F_CHANGESTAT_FN(f_m4_term){
  GET_STORAGE(StoreModelAndStats, storage);

  Free(storage->stats);
  ModelDestroy(nwp, storage->m);
}
