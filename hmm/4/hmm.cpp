#include "hmm.h"

HMM::HMM(int hiddenStates, int emissions) {
  transitionMatrix = vector<vector<long double>>(hiddenStates, vector<long double>(hiddenStates));
  emissionMatrix = vector<vector<long double>>(hiddenStates, vector<long double>(emissions));
  initialMatrix = vector<vector<long double>>(1, vector<long double>(hiddenStates));
  
  // default_random_engine generator((unsigned int)time(0));
  // normal_distribution<long double> distribution(0, 1);

  for (int i = 0; i < hiddenStates; ++i) {
    for (int j = 0; j < hiddenStates; ++j) {
      transitionMatrix[i][j] = 1.0 / hiddenStates;
    }
  }

  for (int i = 0; i < hiddenStates; ++i) {
    // long double delta = 1.0 / 9;
    for (int j = 0; j < emissions; ++j) {
      emissionMatrix[i][j] = 1.0 / emissions;
      // emissionMatrix[i][j] = delta;
      // delta -= abs(distribution(generator) * 0.0001);
    }
  }

  for (int i = 0; i < hiddenStates; ++i) {
    initialMatrix[0][i] = 1.0 / hiddenStates;
  }
}

HMM::HMM(vector<vector<long double>> _transitionMatrix, vector<vector<long
    double>> _emissionMatrix, vector<vector<long double>> _initialMatrix) {
  transitionMatrix = _transitionMatrix;
  emissionMatrix = _emissionMatrix;
  initialMatrix = _initialMatrix;
}

vector<long double> HMM::getNextEmissionDist() {
  // Calculate distribution for all next states.
  vector<long double> nextStateDist(transitionMatrix.size());
  for (unsigned int i = 0; i < transitionMatrix.size(); ++i) {
    for (unsigned int j = 0; j < transitionMatrix[0].size(); ++j) {
      nextStateDist[i] += (initialMatrix[0][j] * transitionMatrix[j][i]);
    }
  }

  // Calculate distribution for all next emissions.
  vector<long double>nextEmissionDist(emissionMatrix[0].size());
  for (unsigned int i = 0; i < emissionMatrix[0].size(); ++i) {
    for (unsigned int j = 0; j < emissionMatrix.size(); ++j) {
      nextEmissionDist[i] += nextStateDist[j] * emissionMatrix[j][i];
    }
  }
  return nextEmissionDist;
}

vector<vector<long double>> HMM::getAlpha(vector<int> &sequence) {
  vector<vector<long double>> alpha(sequence.size(), vector<long double>(transitionMatrix.size()));
  // Initialize. Calculate alphaduder sequence[0].
  for (unsigned int i = 0; i < transitionMatrix.size(); ++i) {
    alpha[0][i] = initialMatrix[0][i] * emissionMatrix[i][sequence[0]];
  }

  // Recurse.
  for (unsigned int t = 1; t < sequence.size(); ++t) {
    for (unsigned int j = 0; j < transitionMatrix.size(); ++j) {
      long double sum = 0;
      for (unsigned int k = 0; k < transitionMatrix[0].size(); ++k) {
        sum += alpha[t - 1][k] * transitionMatrix[k][j];
      }
      alpha[t][j] = sum * emissionMatrix[j][sequence[t]];
    }
  }

  return alpha;
}

vector<vector<long double>> HMM::getBeta(vector<int> &sequence) {
  vector<vector<long double>> beta(sequence.size(), vector<long double>(transitionMatrix.size()));
  // Init.
  for (unsigned int i = 0; i < transitionMatrix.size(); ++i) {
    beta[sequence.size() - 1][i] = 1;
  }

  // Recursion.
  for (unsigned int t = sequence.size() - 1; t > 0; --t) {
    for (unsigned int i = 0; i < transitionMatrix.size(); ++i) {
      long double sum = 0;
      for (unsigned int j = 0; j < transitionMatrix.size(); ++j) {
        sum += transitionMatrix[i][j] * emissionMatrix[j][sequence[t]] * beta[t][j];
      }
      beta[t - 1][i] = sum;
    }
  }
  return beta;
}

long double HMM::getNorm(vector<vector<long double>> &alpha, vector<vector<long double>> &beta) {
  long double sum = 0;
  for (unsigned int i = 0; i < alpha[0].size(); ++i) {
    // cout << "alpha0i: " << alpha[0][i] << " beta0i: " << beta[0][i] << endl;
    sum += alpha[0][i] * beta[0][i];
  }
  return sum;
}

vector<vector<long double>> HMM::getGamma(vector<vector<long double>> &alpha, vector<vector<long double>> &beta, long double norm) {
  vector<vector<long double>> gamma(alpha.size(), vector<long double>(alpha[0].size()));
  for (unsigned int t = 0; t < alpha.size(); ++t) {
    // long double tsum = 0;
    for (unsigned int i = 0; i < alpha[0].size(); ++i) {
      gamma[t][i] = (alpha[t][i]*beta[t][i]) / norm;
      // cout << gamma[t][i] << endl;
      // assert(gamma[t][i] <= 1);
      // assert(gamma[t][i] >= 0);
      // tsum += gamma[t][i];
    }
    // cout << "tsum: " << tsum << endl;
    // assert(abs(tsum - 1) < FAULT);
  }
  return gamma;
}

vector<vector<long double>> HMM::getInitials(vector<vector<long double>> &gamma) {
  vector<vector<long double>> initials(1, vector<long double>(gamma.size()));
  for (unsigned int i = 0; i < gamma.size(); ++i) {
    initials[0][i] = gamma[0][i];
  }
  return initials;
}

long double HMM::getForwardProbability(vector<int> &sequence) {
  vector<vector<long double>> alpha = getAlpha(sequence);

  // Termination.
  long double sum = 0;
  for (unsigned int i = 0; i < transitionMatrix.size(); ++i) {
    sum += alpha[sequence.size() - 1][i];
  }
  return sum;
}

vector<int> HMM::getLikeliestHiddenStates(vector<int> &sequence) {
  vector<vector<long double>> delta(sequence.size(), vector<long double>(transitionMatrix.size()));
  vector<vector<int>> psi(sequence.size(), vector<int>(transitionMatrix.size()));
  vector<int> likeliestHiddenStates(sequence.size());

  // Initialize.
  for (unsigned int i = 0; i < transitionMatrix.size(); ++i) {
    delta[0][i] = initialMatrix[0][i] * emissionMatrix[i][sequence[0]];
    psi[0][i] = 0;
  }

  // Recurse.
  for (unsigned int i = 1; i < sequence.size(); ++i) {
    for (unsigned int j = 0; j < transitionMatrix.size(); ++j) {
      long double p = 0;
      int m = 0;
      for (unsigned int k = 0; k < transitionMatrix.size(); ++k) {
        long double q = delta[i - 1][k] * transitionMatrix[k][j];
        if (q > p) {
          p = q;
          m = k;
        }
      }
      delta[i][j] = p * emissionMatrix[j][sequence[i]];
      psi[i][j] = m;
    }
  }

  // Terminate.
  long double p = 0;
  int m = 0;
  for (unsigned int i = 0; i < transitionMatrix.size(); ++i) {
    long double q = delta[sequence.size() - 1][i];
    if (q > p) {
      p = q;
      m = i;
    }
  }
  likeliestHiddenStates[sequence.size() - 1] = m;

  // Backtrack like void.
  for (unsigned int i = sequence.size() - 1; i > 0; --i) {
    likeliestHiddenStates[i - 1] = psi[i][likeliestHiddenStates[i]];
  }

  return likeliestHiddenStates;
}

void HMM::estimateMatrices(vector<int> &sequence) {
  vector<vector<long double>> oldTransition(transitionMatrix);
  vector<vector<long double>> oldEmission(emissionMatrix);
  long double oldDelta = 100;
  for (int x = 0; x < 25; ++x) {
    // printEmissionMatrix();
    baumWelchIteration(sequence);
    long double deltaSum = 0;
    for (unsigned int i = 0; i < transitionMatrix.size(); ++i) {
      for (unsigned int j = 0; j < transitionMatrix.size(); ++j) 
        deltaSum += abs(transitionMatrix[i][j] - oldTransition[i][j]);
      for (unsigned int j = 0; j < emissionMatrix[0].size(); ++j) 
        deltaSum += abs(emissionMatrix[i][j] - oldEmission[i][j]);
    }
    // cerr << abs(oldDelta - deltaSum) << std::endl;
    if (abs(oldDelta - deltaSum) < 1e-18) {
      // cerr << "BREAKING AT " << x << std::endl;
      break;
    }
    oldDelta = deltaSum;
  }
}

void HMM::baumWelchIteration(vector<int> &sequence) {
  vector<vector<long double>> alpha = getAlpha(sequence);
  vector<vector<long double>> beta = getBeta(sequence);
  long double norm = getNorm(alpha, beta);
  vector<vector<long double>> gamma = getGamma(alpha, beta, norm);
  vector<vector<long double>> newInitialMatrix = getInitials(gamma);


  // Improve transition matrix.
  for (unsigned int i = 0; i < transitionMatrix.size(); ++i) {
    for (unsigned int j = 0; j < transitionMatrix[0].size(); ++j) {
      long double num, den;
      den = num = 0;
      for (unsigned int t = 0; t < sequence.size() - 1; ++t) {
        long double xi = (alpha[t][i] * transitionMatrix[i][j] * emissionMatrix[j][sequence[t + 1]] * beta[t + 1][j]) / norm;
        num += xi;
        den += gamma[t][i];
      }
      // newTransitionMatrix[i][j] = num / den;
      transitionMatrix[i][j] = num / den;
    }
  }

  // Improve emission matrix.
  for (unsigned int i = 0; i < emissionMatrix.size(); ++i) {
    for (unsigned int j = 0; j < emissionMatrix[0].size(); ++j) {
      long double den, num;
      den = num = 0;
      for (unsigned int t = 0; t < sequence.size(); ++t) {
        if ((unsigned int)sequence[t] == j) 
          num += gamma[t][i];
        den += gamma[t][i];
      }
      emissionMatrix[i][j] = num / den;
    }
  }
}

void HMM::printHMM1() {
  std::vector<long double> nextEmissionDist = getNextEmissionDist();

  // Print that shit. HMM1.
  std::cout << 1 << " " << emissionMatrix[0].size();
  for (unsigned int i = 0; i < emissionMatrix[0].size(); ++i) {
    std::cout << " " << nextEmissionDist[i];
  }
}

void HMM::printHMM4() {
  int precision = 6;
  // Print the shit.
  cout << setprecision(precision) <<  transitionMatrix.size() << " " << transitionMatrix[0].size();
  for (unsigned int i = 0; i < transitionMatrix.size(); ++i) {
    for (unsigned int j = 0; j < transitionMatrix[0].size(); ++j) {
      cout << " " << setprecision(precision) << transitionMatrix[i][j];
    }
  }

  cout << endl;

  cout << setprecision(precision) << emissionMatrix.size() << " " << emissionMatrix[0].size();
  for (unsigned int i = 0; i < emissionMatrix.size(); ++i) {
    for (unsigned int j = 0; j < emissionMatrix[0].size(); ++j) {
      cout << " " << setprecision(precision) << emissionMatrix[i][j];
    }
  }
}

void HMM::setTransitionMatrix(vector<vector<long double>> _transitionMatrix) {
  transitionMatrix = _transitionMatrix;
}

void HMM::setEmissionMatrix(vector<vector<long double>> _emissionMatrix) {
  emissionMatrix = _emissionMatrix;
}

void HMM::setInitialMatrix(vector<vector<long double>> _initialMatrix) {
  initialMatrix = _initialMatrix;
}

void HMM::printEmissionMatrix() {
  for (unsigned int i = 0; i < emissionMatrix.size(); ++i) {
    for (unsigned int j = 0; j < emissionMatrix[0].size(); ++j) {
      cerr << emissionMatrix[i][j] << " ";
    }
    cerr << endl;
  }
  cerr << endl;
}
