typedef struct graph{
  int size;

  double phi;
  double eta;
  double tau;
  double lambda;
 
  double* nodeCmf;
  
  int infected;
  int susceptible;
  double avginfected; 

  int **adjmat;
  int *nodevec;
  int *infvec;
}graph;

graph* createGraph( int size, int m0, double phi, double eta);
void setInfVec( graph* g );
void setStatus( graph* g, int node, int status);
void addEdge( graph* g, int nodeA, int nodeB);
void removeEdge( graph* g, int nodeA, int nodeB);
void removeNode( graph* g, int node );
void printGraph( graph* g);
void simEpidemic( graph* g, int burnin, int timesteps );
