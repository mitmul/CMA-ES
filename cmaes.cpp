#include "cmaes.h"

CmaEs::CmaEs(const int _dimension)
  : dimension(_dimension),
    gen(static_cast<unsigned long>(time(0))),
    nd(0.0, 1.0),
    mt(gen, nd)
{
  // 初期平均ベクトル
  xmean.resize(dimension);
  for(int i = 0; i < dimension; ++i)
  {
    xmean[i] = mt();//(double)rand() / RAND_MAX * 2.0 - 1.0;
  }

  // 初期分散
  sigma = 0.5;

  // 収束条件
  stopfitness = 1E-10;
  stopeval = 1E3 * pow(dimension, 2.0);

  // 選択パラメータ
  lambda = 4 + (300 * log(dimension)); // 子孫の数
  mu = lambda / 2;

  // 重みベクトル
  weights.resize(mu);
  weights.setZero();
  for(int i = 0; i < mu; ++i)
  {
    weights[i] = log(mu + 1.0 / 2.0) - log(i + 1.0);
  }
  double weights_sum = weights.sum();
  weights /= weights_sum;
  weights_sum = weights.sum();

  double weights_sum2 = 0.0;
  for(int i = 0; i < mu; ++i)
  {
    weights_sum2 += pow(weights[i], 2.0);
  }
  mueff = pow(weights_sum, 2.0) / weights_sum2;

  // 適応パラメータ
  cc = (4 + mueff / dimension) / (dimension + 4 + 2 * mueff / dimension);
  cs = (mueff + 2.0) / (dimension + mueff + 5.0);
  c1 = 2 / (pow(dimension + 1.3, 2) + mueff);
  cmu = 2 * (mueff - 2 + 1 / mueff) / (pow(dimension + 2, 2.0) + 2 * mueff / 2);
  damps = 1.0 + 2.0 * max(0.0, sqrt((mueff - 1.0) / (dimension + 1.0)) - 1.0) + cs;

  // 動的パラメータ
  pc.resize(dimension);
  pc.setZero();
  ps.resize(dimension);
  ps.setZero();

  // 固有ベクトル行列
  B.resize(dimension, dimension);
  B = MatrixXd::Identity(dimension, dimension);

  // 固有値ベクトル
  D.resize(dimension);
  for(int i = 0; i < dimension; ++i)
  {
    D[i] = 1;
  }

  // 共分散行列
  VectorXd D2 = D;
  for(int i = 0; i < dimension; ++i)
  {
    D2[i] = pow(D[i], 2.0);
  }
  MatrixXd DD(dimension, dimension);
  DD.setZero();
  DD.diagonal() = D2;
  C = B * DD * B.transpose();

  for(int i = 0; i < dimension; ++i)
  {
    D2[i] = pow(D[i], -1.0);
  }
  DD.setZero();
  DD.diagonal() = D2;
  invsqrtC = B * DD * B.transpose();

  eigeneval = 0;
  chiN = pow(dimension, 0.5) * (1.0 - 1.0 / (4.0 * dimension) + 1.0 / (21.0 * pow(dimension, 2.0)));

  arx.resize(dimension, lambda);
  arfitness.resize(lambda);
  counteval = 0;
  arindex.resize(lambda);
}

CmaEs::~CmaEs()
{
}

double CmaEs::cost(const VectorXd &input_x)
{
#define ROSENBROCK

#ifdef RASTRIGIN
  int n = input_x.rows();

  VectorXd scale(n);
  for(int i = 0; i < n; ++i)
    scale[i] = pow(10.0, (double)i / ((double)n - 1.0));

  double sum = 0.0;
  for(int i = 0; i < n; ++i)
    sum += pow(scale[i] * input_x[i], 2.0) - 10.0 * cos(2.0 * M_PI * scale[i] * input_x[i]);

  return 10.0 * n + sum;
#endif

#ifdef ROSENBROCK
  double sum = 0.0;
  int n = input_x.rows();
  for(int i = 0; i < n - 1; ++i)
  {
    sum += 100.0 * pow(input_x[i + 1] - pow(input_x[i], 2.0), 2.0) + pow(1.0 - input_x[i], 2.0);
  }

  return sum;
#endif
}

void CmaEs::generateOffspring()
{
  arx.setZero();

  for(int i = 0; i < lambda; ++i)
  {
    VectorXd D_rand = D;
    for(int j = 0; j < dimension; ++j)
    {
      D_rand[j] *= mt();
    }
    arx.col(i) = xmean + sigma * B * D_rand;

    arfitness[i] = cost(arx.col(i));
    ++counteval;
  }

  map<double, int> arfitness_map;
  for(int i = 0; i < lambda; ++i)
  {
    arfitness_map.insert(pair<double, int>(arfitness[i], i));
  }

  int i = 0;
  arindex.resize(lambda);
  for(map<double, int>::iterator it = arfitness_map.begin(); it != arfitness_map.end(); ++it, ++i)
  {
    arfitness[i] = (*it).first;
    arindex[i] = (*it).second;
  }
}

void CmaEs::updateMean()
{
  MatrixXd good_arx(dimension, mu);
  for(int i = 0; i < mu; ++i)
    good_arx.col(i) = arx.col(arindex[i]);

  xold = xmean;
  xmean = good_arx * weights;
}

void CmaEs::updatePath()
{
  ps = (1.0 - cs) * ps + sqrt(cs * (2.0 - cs) * mueff) * invsqrtC * (xmean - xold) / sigma;
  hsig = ps.squaredNorm() / (1.0 - pow(1.0 - cs, 2.0 * (double)counteval / (double)lambda)) / dimension < (2.0 + 4.0 / (dimension + 1.0)) ? 1 : 0;
  pc = (1.0 - cc) * pc + hsig * sqrt(cc * (2.0 - cc) * mueff) * (xmean - xold) / sigma;
}

void CmaEs::adaptCovmat()
{
  MatrixXd good_arx(dimension, mu);
  for(int i = 0; i < mu; ++i)
    good_arx.col(i) = arx.col(arindex[i]);

  MatrixXd xold_rep(dimension, mu);
  for(int i = 0; i < mu; ++i)
    xold_rep.col(i) = xold;

  MatrixXd artmp(dimension, mu);
  artmp = (1.0 / sigma) * (good_arx - xold_rep);

  MatrixXd weights_diag(mu, mu);
  weights_diag.setZero();
  weights_diag.diagonal() = weights;

  C = (1.0 - c1 - cmu) * C + c1 * (pc * pc.transpose() + (1 - hsig) * cc * (2.0 - cc) * C) + cmu * artmp * weights_diag * artmp.transpose();
}

void CmaEs::adaptSigma()
{
  sigma = sigma * exp((cs / damps) * (ps.norm() / chiN - 1.0));
}

void CmaEs::updateEigen()
{
  if(counteval - eigeneval > lambda / (c1 + cmu) / dimension / 10.0)
  {
    eigeneval = counteval;

    MatrixXd triu_C = C;
    for(int i = 0; i < dimension; ++i)
      for(int j = i + 1; j < dimension; ++j)
        triu_C(j, i) = 0;

    MatrixXd triu_C1 = C;
    for(int i = 0; i < dimension; ++i)
      for(int j = i; j < dimension; ++j)
        triu_C1(j, i) = 0;

    C = triu_C + triu_C1.transpose();

    SelfAdjointEigenSolver<MatrixXd> eigensolver(C);
    VectorXd val = eigensolver.eigenvalues();
    for(int i = 0; i < dimension; ++i)
      val[i] = sqrt(val[i]);
    D = val;
    B = eigensolver.eigenvectors();

    MatrixXd D_diag(dimension, dimension);
    D_diag.setZero();
    for(int i = 0; i < dimension; ++i)
      D_diag(i, i) = 1.0 / D[i];

    invsqrtC = B * D_diag * B.transpose();
  }
}

void CmaEs::loadMatrix(MatrixXd &mat, const string &file_name)
{
  vector<VectorXd> vecs;
  ifstream ifs(file_name.c_str());
  string buffer;
  while(ifs && getline(ifs, buffer))
  {
    list<string> results;
    boost::split(results, buffer, boost::is_any_of(" "));

    int i = 0;
    VectorXd vec(results.size());
    for(list<string>::iterator it = results.begin(); it != results.end(); ++it)
    {
      string st = (*it);
      if(st == "#")
        break;
      else if(st == " " || st.empty())
        continue;
      vec[i] = boost::lexical_cast<double>(st);
      ++i;
    }
    if(i != 0)
      vecs.push_back(vec);
  }

  mat.resize(vecs.size(), vecs[0].rows());

  for(int i = 0; i < (int)vecs.size(); ++i)
  {
    mat.row(i) = vecs[i];
  }
}

void CmaEs::generationLoop()
{
  counteval = 0;
  while(counteval < stopeval)
  {
    generateOffspring();
    updateMean();
    updatePath();
    adaptCovmat();
    adaptSigma();
    updateEigen();

    double max_d = D[0];
    double min_d = D[0];
    for(int i = 0; i < dimension; ++i)
    {
      if(max_d < D[i])
        max_d = D[i];
      if(min_d > D[i])
        min_d = D[i];
    }

    if(arfitness[0] <= stopfitness || max_d > 1E7 * min_d)
    {
      cout << "arfitness[0]: " << arfitness[0] << endl;
      cout << "stopfitness: " << stopfitness << endl << endl;
      break;
    }

    cout << "xmean:\t";
    for(int i = 0; i < dimension; ++i)
      printf("%.5f\t", xmean[i]);
    cout << endl;

    cout << counteval << ":\t";
    printf("cost:%.5f\t\n", arfitness[0]);
  }
}
