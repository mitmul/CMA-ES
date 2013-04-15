#ifndef CMAES_H
#define CMAES_H

#include <map>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <Eigen/Eigen>
#include <boost/random.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace Eigen;
using namespace boost::random;

class CmaEs
{
  public:
    CmaEs(const int _dimension, const vector<double> _scaling = vector<double>());
    void initialize();

    virtual double cost(const VectorXd &input_x);

    void generationLoop();

  private:
    int dimension;
    vector<double> scaling;

    VectorXd xmean;
    double sigma;
    double stopfitness;
    double stopeval;
    int lambda;
    int mu;
    VectorXd weights;
    double mueff;
    double cc;
    double cs;
    double c1;
    double cmu;
    double damps;
    VectorXd pc;
    VectorXd ps;
    MatrixXd B;
    VectorXd D;
    MatrixXd C;
    MatrixXd invsqrtC;
    double eigeneval;
    double chiN;

    MatrixXd arx;
    VectorXd arfitness;

    int counteval;
    VectorXi arindex;

    VectorXd xold;

    int hsig;

    void generateOffspring();
    void updateMean();
    void updatePath();
    void adaptCovmat();
    void adaptSigma();
    void updateEigen();

    // 乱数生成用
    mt19937 gen;
    normal_distribution<> nd;
    uniform_real_distribution<> ud;
    variate_generator<mt19937, normal_distribution<> > randn;
    variate_generator<mt19937, uniform_real_distribution<> > randu;

    void loadMatrix(MatrixXd &mat, const string &file_name);
};

#endif // CMAES_H
