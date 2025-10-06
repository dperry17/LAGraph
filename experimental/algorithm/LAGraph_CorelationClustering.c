#include "LAGraphX.h"
#include "LG_internal.h"

#undef LG_FREE_WORK
#define LG_FREE_WORK                                                           \
  {                                                                            \
    GrB_free(&I);                                                              \
    GrB_free(&pivot_id);                                                       \
    GrB_free(&theta);                                                          \
    GrB_free(&X);                                                              \
    GrB_free(&val);                                                            \
  }

#undef LG_FREE_ALL
#define LG_FREE_ALL				\
  LG_FREE_WORK


int LG_compute_objective(GrB_Matrix A, GrB_Vector C, const GrB_Index n, const double plus_weight, const double minus_weight, double* score, char* msg) {

  GrB_Vector I = NULL, pivot_id = NULL ;
  GrB_Matrix theta = NULL, X = NULL ;
  GrB_Scalar val = NULL ;
  GrB_Index n_clustered = 0, n_unclustered = 0 ;

  GRB_TRY(GrB_Vector_new(&I, GrB_INT64, n)) ;
  GRB_TRY(GrB_Vector_new(&pivot_id, GrB_INT64, n)) ;
  GRB_TRY(GrB_Matrix_new(&theta, GrB_BOOL, n, n)) ;
  GRB_TRY(GrB_Matrix_new(&X, GrB_BOOL, n, n)) ;
  GRB_TRY(GrB_Scalar_new(&val, GrB_BOOL)) ;
  GRB_TRY(GrB_Scalar_setElement(val, true)) ;
  
  GRB_TRY(GxB_Vector_extractTuples_Vector(I, pivot_id, C, NULL)) ;
  GRB_TRY(GrB_Matrix_build(theta, pivot_id, I, val, NULL)) ; //n by n matrix with n values

  //X<A, struct> = A .^ theta, where X holds the unclustered edges
  GRB_TRY(GrB_eWiseMult(X, A, NULL, GxB_LXOR_INT64, A, theta, GrB_DESC_S)) ;

  GRB_TRY(GrB_Matrix_nvals(&n_unclustered, X)) ;
  GRB_TRY(GrB_Matrix_nvals(&n_clustered, theta)) ;

  LG_FREE_ALL;
  *score = (n_clustered * plus_weight) + (n_unclustered * minus_weight) ;

  return (GrB_SUCCESS) ;
}


#undef LG_FREE_WORK
#define LG_FREE_WORK				                                \
{						                                \
 GrB_free(&pivots) ;					                        \
 GrB_free(&clustering) ;						        \
 GrB_free(&create_tuple) ;					                \
 GrB_free(&create_tuple_indexOp) ;						\
 GrB_free(&thunk) ;							        \
 GrB_free(&Seed) ;							        \
 GrB_free(&tuple_min) ;			                                        \
 GrB_free(&tuple_min_monoid) ;			                                \
}

#undef LG_FREE_ALL
#define LG_FREE_ALL                                                            \
{                                                                              \
    LG_FREE_WORK							       \
}

#define F_UNARY(f) ((void (*)(void *, const void *))f)
#define F_BINARY(f) ((void (*)(void *, const void *, const void *))f)
#define F_INDEX_BINARY(f) ((void (*)(void*, const void*, GrB_Index, GrB_Index, const void *, GrB_Index, GrB_Index, const void *)) f)

#define JIT_STR(func, var)                                                     \
  char *var = #func;                                                           \
  func

JIT_STR(
    typedef struct {
      int64_t i;
      double seed;
    } LG_neighbor_tuple;
,NTS)    


JIT_STR(
    void LG_createTuple(LG_neighbor_tuple *z, const uint64_t *x, uint64_t ix,
                        uint64_t jx, const uint64_t *y, uint64_t iy, uint64_t jy,
                        const bool *thunk) {
      
      z->i = iy;
      z->seed = *x;
    },
CC_CTUPLE)

JIT_STR(
	void LG_tupleMin(LG_neighbor_tuple* z, const LG_neighbor_tuple* x, const LG_neighbor_tuple* y){
	  *z = (x->seed < y->seed) ? *x : *y ;
	},
CC_MINT)

JIT_STR(
	void LG_getId(int64_t* z, const LG_neighbor_tuple* x){
	  *z = x->i;
	},
CC_ID)


int LAGraph_CorrelationClustering(GrB_Vector* clusters, const LAGraph_Graph G, char* msg) {

  GrB_Vector pivots = NULL, Seed = NULL, tuple_vector = NULL ;
  GrB_Semiring clustering = NULL ;
  GrB_BinaryOp tuple_min = NULL, create_tuple = NULL ;
  GrB_UnaryOp get_id = NULL ;
  GxB_IndexBinaryOp create_tuple_indexOp = NULL ;
  GrB_Monoid tuple_min_monoid = NULL ;
  GrB_Scalar thunk=NULL, empty = NULL ;
  GrB_Type neighbor_tuple = NULL ;

  GrB_Matrix A = G->A ;

  GrB_Index n = 0;

  GRB_TRY(GxB_Type_new(&neighbor_tuple, sizeof(LG_neighbor_tuple), "LG_neighbor_tuple", NTS)) ;
  GRB_TRY(GrB_Matrix_nrows(&n, A)) ;

  //make graph symmetric
  GRB_TRY(GrB_assign(A, A, NULL, A, GrB_ALL, n, GrB_ALL, n, GrB_DESC_SCT0)) ;

  GRB_TRY(GrB_Vector_new(&pivots, GrB_INT64, n)) ;
  GRB_TRY(GrB_Vector_new(&Seed, GrB_UINT64, n)) ;
  GRB_TRY(GrB_Vector_new(&tuple_vector, neighbor_tuple, n)) ;
  GRB_TRY(GrB_Scalar_new(&empty, GrB_BOOL)) ;

  
  GRB_TRY(GrB_Vector_new(clusters, GrB_INT64, n)) ;
  
  LG_TRY(LAGraph_Cached_OutDegree(G, msg)) ;
  LG_TRY(LAGraph_Cached_NSelfEdges(G, msg)) ;
  LG_TRY(LAGraph_MaximalIndependentSet(&pivots, G, 0, NULL, msg)) ;
    
  //create semiring
  GRB_TRY(GrB_Scalar_new(&thunk, GrB_BOOL)) ;
  GRB_TRY(GrB_Scalar_setElement(thunk, true)) ;
  GRB_TRY(GxB_IndexBinaryOp_new(&create_tuple_indexOp, F_INDEX_BINARY(LG_createTuple),
  				neighbor_tuple, GrB_UINT64, GrB_INT64, GrB_BOOL, "LG_createTuple", CC_CTUPLE)) ;
  GRB_TRY(GxB_BinaryOp_new_IndexOp(&create_tuple, create_tuple_indexOp, thunk)) ;
  GRB_TRY(GxB_BinaryOp_new(&tuple_min, F_BINARY(LG_tupleMin), neighbor_tuple, neighbor_tuple, neighbor_tuple, "LG_tupleMin", CC_MINT)) ;
  GRB_TRY(GxB_UnaryOp_new(&get_id, F_UNARY(LG_tupleMin), GrB_INT64, neighbor_tuple, "LG_getId", CC_ID)) ;
  LG_neighbor_tuple id = {0, 0} ;
  GRB_TRY(GrB_Monoid_new_UDT(&tuple_min_monoid, tuple_min, &id)) ;
  GRB_TRY(GrB_Semiring_new(&clustering, tuple_min_monoid, create_tuple)) ;
  //GxB_print(pivots, 5);
  
  //get all clusters through single mxv operation
  GRB_TRY(GrB_vxm(tuple_vector, pivots, NULL, clustering, pivots, A, GrB_DESC_SC)) ;
  GRB_TRY(GrB_apply(*clusters, pivots, NULL, get_id, tuple_vector, GrB_DESC_SC)) ;
  GRB_TRY (GrB_assign (*clusters, pivots, NULL, 0, GrB_ALL, n, GrB_DESC_S)) ;
  GRB_TRY(GrB_apply(*clusters, pivots, NULL, GrB_ROWINDEX_INT64, *clusters, 0, GrB_DESC_S)) ;

  double score = 0;
  LG_TRY(LG_compute_objective(A, *clusters, n, -1, 1, &score, msg)) ;

  printf("The objective score is: %lf\n", score);
  
  LG_FREE_ALL;
  return (GrB_SUCCESS) ;
  
}
