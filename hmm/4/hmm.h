#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <cmath>
#include <cassert>

#define FAULT 0.0000001

using namespace std;

class HMM {
  private:
    vector<vector<long double>> transitionMatrix;
    vector<vector<long double>> emissionMatrix;
    vector<vector<long double>> initialMatrix;
  public:
    HMM(int hiddenStates, int emissions);
    HMM(vector<vector<long double>> _transitionMatrix, vector<vector<long
        double>> _emissionMatrix, vector<vector<long double>> _initialMatrix);
    vector<long double> getNextEmissionDist();
    vector<vector<long double>> getAlpha(vector<int> &sequence);
    vector<vector<long double>> getBeta(vector<int> &sequence);
    long double getNorm(vector<vector<long double>> &alpha, vector<vector<long double>> &beta);
    vector<vector<long double>> getGamma(vector<vector<long double>> &alpha, vector<vector<long double>> &beta, long double norm);
    vector<vector<long double>> getInitials(vector<vector<long double>> &gamma);
    long double getForwardProbability(vector<int> &sequence);
    vector<int> getLikeliestHiddenStates(vector<int> &sequence);
    void estimateMatrices(vector<int> &sequence);
    void baumWelchIteration(vector<int> &sequence);
    void printHMM1();
    void printHMM4();

    void setTransitionMatrix(vector<vector<long double>> _transitionMatrix);
    void setEmissionMatrix(vector<vector<long double>> _emissionMatrix);
    void setInitialMatrix(vector<vector<long double>> _initialMatrix);
};
