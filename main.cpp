#include <iostream>
#include "cmaes.h"

using namespace std;

int main()
{
  CmaEs *cma_es = new CmaEs(10);
  cma_es->generationLoop();

  cout << "Hello World!" << endl;
  return 0;
}


