#include <iostream>
#include "cmaes.h"

using namespace std;

int main()
{
  int dimension = 20;
  vector<double> scale;
  for(int i = 0; i < dimension; ++i)
  {
    scale.push_back(512);
  }
  CmaEs *cma_es = new CmaEs(dimension, scale);
  cma_es->generationLoop();

  cout << "Hello World!" << endl;
  return 0;
}


