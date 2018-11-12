#include "ergm_changestat_operator.h"

typedef struct{Model *m; double *stats;} StoreModelAndStats;

I_CHANGESTAT_FN(i_taper_term){
  double *inputs = INPUT_PARAM + N_CHANGE_STATS - 1; // Skip the penalty coefficients.
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
      UPDATE_STORAGE(tail, head, nw, m, NULL);
      ToggleEdge(tail, head, nw);
    });
  
  // Note that nw is now identitical to nwp.
  NetworkDestroy(nw);
}

D_CHANGESTAT_FN(d_taper_term){
  GET_STORAGE(StoreModelAndStats, storage);
  Model *m = storage->m;

  // TODO: cache the result for use in u_taper_term.
  ChangeStats(ntoggles, tails, heads, nwp, m);
  
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

  // TODO: Use the result cached by the d_ function.
  ChangeStats(1, &tail, &head, nwp, m);
  for(unsigned int k=0; k<m->n_stats; k++)
    storage->stats[k] += m->workspace[k];

  UPDATE_STORAGE(tail, head, nwp, m, NULL);
}

F_CHANGESTAT_FN(f_taper_term){
  GET_STORAGE(StoreModelAndStats, storage);

  Free(storage->stats);
  ModelDestroy(nwp, storage->m);
}
