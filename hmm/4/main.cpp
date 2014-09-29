#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <cmath>

std::vector<long double> getNextEmissionDist(std::vector<std::vector<long double>> transitionMatrix, 
                                        std::vector<std::vector<long double>> emissionMatrix, 
                                        std::vector<std::vector<long double>> initialMatrix);
std::vector<std::vector<long double>> getAlpha(std::vector<std::vector<long double>> transitionMatrix, 
                                          std::vector<std::vector<long double>> emissionMatrix, 
                                          std::vector<std::vector<long double>> initialMatrix,
                                          std::vector<int> sequence);
std::vector<std::vector<long double>> getBeta(std::vector<std::vector<long double>> transitionMatrix, 
                                          std::vector<std::vector<long double>> emissionMatrix, 
                                          std::vector<std::vector<long double>> initialMatrix,
                                          std::vector<int> sequence);
long double getNorm(std::vector<std::vector<long double>> alpha, std::vector<std::vector<long double>> beta);
std::vector<std::vector<long double>> getGamma(std::vector<std::vector<long double>> alpha, std::vector<std::vector<long double>> beta, long double norm);
long double getForwardProbability(std::vector<std::vector<long double>> transitionMatrix,
                                        std::vector<std::vector<long double>> emissionMatrix,
                                        std::vector<std::vector<long double>> initialMatrix,
                                        std::vector<int> sequence);
std::vector<int> getLikeliestHiddenStates(std::vector<std::vector<long double>> transitionMatrix,
                                          std::vector<std::vector<long double>> emissionMatrix,
                                          std::vector<std::vector<long double>> initialMatrix,
                                          std::vector<int> sequence);

// Most of the stuff in this file is stolen from Nikolai Shokhirev.
// His code was written in Pascal though, so some things might differ.
// http://www.shokhirev.com/nikolai/abc/alg/hmm/hmm.html
int main(int argc, char *argv[]) {
  // Read stuff from cin. 
  std::string transitionString;
  std::getline(std::cin, transitionString);
  
  std::string emissionString;
  std::getline(std::cin, emissionString);

  std::string initialString;
  std::getline(std::cin, initialString);

  std::string sequenceString;
  std::getline(std::cin, sequenceString);

  // Transitions.
  int transitionRows, transitionColumns;
  std::istringstream iss(transitionString);
  iss >> transitionRows;
  iss >> transitionColumns;

  int x, y; 
  long double tmp;
  std::vector<std::vector<long double>> transitions(transitionRows, std::vector<long double>(transitionColumns));
  x = y = 0;
  while (iss >> tmp) {
    transitions[x][y++] = tmp;
    if (y == transitionColumns) { ++x; y = 0; }
  }

  // Emissions.
  int emissionRows, emissionColumns;
  std::istringstream iss1(emissionString);
  iss1 >> emissionRows;
  iss1 >> emissionColumns;

  std::vector<std::vector<long double>> emissions(emissionRows, std::vector<long double>(emissionColumns));
  x = y = 0;
  while (iss1 >> tmp) {
    emissions[x][y++] = tmp;
    if (y == emissionColumns) { ++x; y = 0; }
  }

  // Initials.
  int initialRows, initialColumns;
  std::istringstream iss2(initialString);
  iss2 >> initialRows;
  iss2 >> initialColumns;

  std::vector<std::vector<long double>> initials(initialRows, std::vector<long double>(initialColumns));
  x = y = 0;
  while (iss2 >> tmp) {
    initials[x][y++] = tmp;
    if (y == initialColumns) { ++x; y = 0; }
  }

  // Sequence.
  int sequenceLength;
  std::istringstream iss3(sequenceString);
  iss3 >> sequenceLength;
  std::vector<int> sequence(sequenceLength);
  x = 0;
  while (iss3 >> tmp) {
    sequence[x++] = tmp;
  }

  std::vector<std::vector<long double>> alpha = getAlpha(transitions, emissions, initials, sequence);
  std::vector<std::vector<long double>> beta = getBeta(transitions, emissions, initials, sequence);
  long double norm = getNorm(alpha, beta);
  std::vector<std::vector<long double>> gamma = getGamma(alpha, beta, norm);

  // Improve transition matrix.
  std::vector<std::vector<long double>> newTransitions(transitions.size(), std::vector<long double>(transitions[0].size()));
  for (unsigned int i = 0; i < transitions.size(); ++i) {
    long double den = 0;
    for (unsigned int t = 0; t < sequence.size() - 1; ++t)
      den += gamma[t][i];

    for (unsigned int j = 0; j < transitions[0].size(); ++j) {
      long double num = 0;
      for (unsigned int t = 0; t < sequence.size() - 1; ++t) {
        long double xi = (alpha[t][i] * transitions[i][j] * emissions[j][sequence[t + 1]] * beta[t + 1][i]) / norm;
        num += xi;
      }
      newTransitions[i][j] = num / den;
    }
  }

  // Improve emission matrix.
  std::vector<std::vector<long double>> newEmissions(emissions.size(), std::vector<long double>(emissions[0].size()));
  for (unsigned int i = 0; i < emissions.size(); ++i) {
    for (unsigned int j = 0; j < emissions[0].size(); ++j) {
      newEmissions[i][j] = 0;
    }
  }

  for (unsigned int i = 0; i < transitions.size(); ++i) {
    long double den = 0;
    for (unsigned int t = 0; t < sequence.size(); ++t) {
      den += gamma[t][i];
      for (unsigned int t = 0; t < sequence.size(); ++t) {
        newEmissions[i][sequence[t]] = (newEmissions[i][sequence[t]] + gamma[t][i]) / den;
      }
    }
  }

  // Print the shit.
  std::cout << newTransitions.size() << " " << newTransitions[0].size();
  for (unsigned int i = 0; i < transitions.size(); ++i) {
    for (unsigned int j = 0; j < transitions[0].size(); ++j) {
      std::cout << " " << newTransitions[i][j];
    }
  }

  std::cout << std::endl;

  std::cout << newEmissions.size() << " " << newEmissions[0].size();
  for (unsigned int i = 0; i < emissions.size(); ++i) {
    for (unsigned int j = 0; j < emissions[0].size(); ++j) {
      std::cout << " " << newEmissions[i][j];
    }
  }
  return 0;
}

std::vector<long double> getNextEmissionDist(std::vector<std::vector<long double>> transitionMatrix, 
                                        std::vector<std::vector<long double>> emissionMatrix, 
                                        std::vector<std::vector<long double>> initialMatrix) {
  // Calculate distribution for all next states.
  std::vector<long double> nextStateDist(transitionMatrix.size());
  for (unsigned int i = 0; i < transitionMatrix.size(); ++i) {
    for (unsigned int j = 0; j < transitionMatrix[0].size(); ++j) {
      nextStateDist[i] += (initialMatrix[0][j] * transitionMatrix[j][i]);
    }
  }

  // Calculate distribution for all next emissions.
  std::vector<long double>nextEmissionDist(emissionMatrix[0].size());
  for (unsigned int i = 0; i < emissionMatrix[0].size(); ++i) {
    for (unsigned int j = 0; j < emissionMatrix.size(); ++j) {
      nextEmissionDist[i] += nextStateDist[j] * emissionMatrix[j][i];
    }
  }
  return nextEmissionDist;
}

std::vector<std::vector<long double>> getAlpha(std::vector<std::vector<long double>> transitionMatrix, 
                                          std::vector<std::vector<long double>> emissionMatrix, 
                                          std::vector<std::vector<long double>> initialMatrix,
                                          std::vector<int> sequence) {
  std::vector<std::vector<long double>> alpha(sequence.size(), std::vector<long double>(transitionMatrix.size()));
  // Initialize. Calculate alphaduder sequence[0].
  for (unsigned int i = 0; i < transitionMatrix.size(); ++i) {
    alpha[0][i] = initialMatrix[0][i] * emissionMatrix[i][sequence[0]];
  }

  // Recurse.
  for (unsigned int i = 1; i < sequence.size(); ++i) {
    for (unsigned int j = 0; j < transitionMatrix.size(); ++j) {
      long double sum = 0;
      for (unsigned int k = 0; k < transitionMatrix[0].size(); ++k) {
        sum += alpha[i - 1][k] * transitionMatrix[k][j];
      }
      alpha[i][j] = sum * emissionMatrix[j][sequence[i]];
    }
  }
  return alpha;
}

std::vector<std::vector<long double>> getBeta(std::vector<std::vector<long double>> transitionMatrix, 
                                          std::vector<std::vector<long double>> emissionMatrix, 
                                          std::vector<std::vector<long double>> initialMatrix,
                                          std::vector<int> sequence) {
  std::vector<std::vector<long double>> beta(sequence.size(), std::vector<long double>(transitionMatrix.size()));
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

long double getNorm(std::vector<std::vector<long double>> alpha, std::vector<std::vector<long double>> beta) {
  long double sum = 0;
  for (unsigned int i = 0; i < alpha[0].size(); ++i) {
    // std::cout << "alpha0i: " << alpha[0][i] << " beta0i: " << beta[0][i] << std::endl;
    sum += alpha[0][i] * beta[0][i];
  }
  return sum;
}

std::vector<std::vector<long double>> getGamma(std::vector<std::vector<long double>> alpha, std::vector<std::vector<long double>> beta, long double norm) {
  std::vector<std::vector<long double>> gamma(alpha.size(), std::vector<long double>(alpha[0].size()));
  for (unsigned int t = 0; t < alpha.size(); ++t) {
    for (unsigned int i = 0; i < alpha[0].size(); ++i) {
      gamma[t][i] = (alpha[t][i]*beta[t][i]) / norm;
    }
  }
  return gamma;
}

long double getForwardProbability(std::vector<std::vector<long double>> transitionMatrix,
                             std::vector<std::vector<long double>> emissionMatrix,
                             std::vector<std::vector<long double>> initialMatrix,
                             std::vector<int> sequence) {
  std::vector<std::vector<long double>> alpha = getAlpha(transitionMatrix, emissionMatrix, initialMatrix, sequence);

  // Termination.
  long double sum = 0;
  for (unsigned int i = 0; i < transitionMatrix.size(); ++i) {
    sum += alpha[sequence.size() - 1][i];
  }
  return sum;
}

std::vector<int> getLikeliestHiddenStates(std::vector<std::vector<long double>> transitionMatrix,
                                          std::vector<std::vector<long double>> emissionMatrix,
                                          std::vector<std::vector<long double>> initialMatrix,
                                          std::vector<int> sequence) {
  std::vector<std::vector<long double>> delta(sequence.size(), std::vector<long double>(transitionMatrix.size()));
  std::vector<std::vector<int>> psi(sequence.size(), std::vector<int>(transitionMatrix.size()));
  std::vector<int> likeliestHiddenStates(sequence.size());

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
