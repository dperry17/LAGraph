#include <acutest.h>
#include <LAGraphX.h>
#include <LAGraph_test.h>
#include <stdio.h>
#include <LG_Xtest.h>
#include "LG_internal.h"


char msg[LAGRAPH_MSG_LEN];
LAGraph_Graph G = NULL;
GrB_Matrix A = NULL;
#define LEN 512
#define NTESTS 5
char filename[LEN + 1];

typedef struct{
  char* filename;
  LAGraph_Kind kind ;
}test_info;

test_info tests[] = {
  {"mcl.mtx", LAGraph_ADJACENCY_UNDIRECTED},
  {"belgium_osm.mtx", LAGraph_ADJACENCY_UNDIRECTED},
  {"rgg_n_2_19_s0.mtx", LAGraph_ADJACENCY_UNDIRECTED},
  {"delaunay_n20.mtx", LAGraph_ADJACENCY_UNDIRECTED},
  {"road_central.mtx", LAGraph_ADJACENCY_UNDIRECTED}
};

void test_CorrelationClustering(void) {
#if LG_SUITESPARSE_GRAPHBLAS_V10
  LAGraph_Init(msg);
//OK(LG_SET_BURBLE(1));
  OK(LG_SET_BURBLE(0));
  int inner_threads = 0, outer_threads = 0;
  LAGraph_GetNumThreads(&outer_threads, &inner_threads, msg) ;
  printf("num inner threads: %d\n", inner_threads) ;
  for(uint8_t test = 0; test < NTESTS; test++){
    GrB_Matrix A=NULL;
    GrB_Vector clusters = NULL ;
    printf ("\nMatrix: %s\n", tests[test].filename);
    TEST_CASE(tests[test].filename);
    snprintf(filename, LEN, LG_DATA_DIR "%s", tests[test].filename);
    FILE* f = fopen(filename, "r");
    
    OK(LAGraph_MMRead(&A, f, msg));
    OK(fclose(f));
    LAGraph_Kind kind = tests [test].kind ;
    OK(LAGraph_New(&G, &A, kind, msg));
   
    // test with JIT
    OK(GxB_Global_Option_set(GxB_JIT_C_CONTROL, GxB_JIT_ON));
    double time = LAGraph_WallClockTime() ;
    OK(LAGraph_CorrelationClustering(&clusters, G, msg));
    time = LAGraph_WallClockTime() - time ;
    printf("%s\n", msg);
    printf("Time for Correlation clustering is: %lf", time) ;
    //GxB_print(clusters, 5);

   
    //free work
    OK(LAGraph_Delete(&G, msg));
    GrB_free(&clusters); 
  }
  LAGraph_Finalize(msg);
#endif
}

TEST_LIST = {{"CorrelationClustering", test_CorrelationClustering}, {NULL, NULL}};
