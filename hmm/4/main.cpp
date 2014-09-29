#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <cmath>

std::vector<double> getNextEmissionDist(std::vector<std::vector<double>> transitionMatrix, 
                                        std::vector<std::vector<double>> emissionMatrix, 
                                        std::vector<std::vector<double>> initialMatrix);
double getProbabilityOfEmissionSequence(std::vector<std::vector<double>> transitionMatrix,
                                        std::vector<std::vector<double>> emissionMatrix,
                                        std::vector<std::vector<double>> initialMatrix,
                                        std::vector<int> sequence);
std::vector<int> getLikeliestHiddenStates(std::vector<std::vector<double>> transitionMatrix,
                                          std::vector<std::vector<double>> emissionMatrix,
                                          std::vector<std::vector<double>> initialMatrix,
                                          std::vector<int> sequence);

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
  double tmp;
  std::vector<std::vector<double>> transitions(transitionRows, std::vector<double>(transitionColumns));
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

  std::vector<std::vector<double>> emissions(emissionRows, std::vector<double>(emissionColumns));
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

  std::vector<std::vector<double>> initials(initialRows, std::vector<double>(initialColumns));
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


  std::vector<int> likeliestHiddenStates = getLikeliestHiddenStates(transitions, emissions, initials, sequence);
  std::cout << likeliestHiddenStates[0];
  for (unsigned int i = 1; i < likeliestHiddenStates.size(); ++i) {
  std::cout << " " << likeliestHiddenStates[i];
  }

  return 0;
}

std::vector<double> getNextEmissionDist(std::vector<std::vector<double>> transitionMatrix, 
                                        std::vector<std::vector<double>> emissionMatrix, 
                                        std::vector<std::vector<double>> initialMatrix) {
  // Calculate distribution for all next states.
  std::vector<double> nextStateDist(transitionMatrix.size());
  for (unsigned int i = 0; i < transitionMatrix.size(); ++i) {
    for (unsigned int j = 0; j < transitionMatrix[0].size(); ++j) {
      nextStateDist[i] += (initialMatrix[0][j] * transitionMatrix[j][i]);
    }
  }

  // Calculate distribution for all next emissions.
  std::vector<double>nextEmissionDist(emissionMatrix[0].size());
  for (unsigned int i = 0; i < emissionMatrix[0].size(); ++i) {
    for (unsigned int j = 0; j < emissionMatrix.size(); ++j) {
      nextEmissionDist[i] += nextStateDist[j] * emissionMatrix[j][i];
    }
  }
  return nextEmissionDist;
}

double getProbabilityOfEmissionSequence(std::vector<std::vector<double>> transitionMatrix,
                                        std::vector<std::vector<double>> emissionMatrix,
                                        std::vector<std::vector<double>> initialMatrix,
                                        std::vector<int> sequence) {
  std::vector<std::vector<double>> alpha(sequence.size(), std::vector<double>(transitionMatrix.size()));
  // Initialize. Calculate alphaduder sequence[0].
  for (unsigned int i = 0; i < transitionMatrix.size(); ++i) {
    alpha[0][i] = initialMatrix[0][i] * emissionMatrix[i][sequence[0]];
  }

  // Recurse.
  for (unsigned int i = 1; i < sequence.size(); ++i) {
    for (unsigned int j = 0; j < transitionMatrix.size(); ++j) {
      double sum = 0;
      for (unsigned int k = 0; k < transitionMatrix.size(); ++k) {
        sum += alpha[i - 1][k] * transitionMatrix[k][j];
      }
      alpha[i][j] = sum * emissionMatrix[j][sequence[i]];
    }
  }

  // Termination.
  double sum = 0;
  for (unsigned int i = 0; i < transitionMatrix.size(); ++i) {
    sum += alpha[sequence.size() - 1][i];
  }
  return sum;
}

// Stolen from Nikolai Shokhirev. http://www.shokhirev.com/nikolai/abc/alg/hmm/hmm.html
std::vector<int> getLikeliestHiddenStates(std::vector<std::vector<double>> transitionMatrix,
                                          std::vector<std::vector<double>> emissionMatrix,
                                          std::vector<std::vector<double>> initialMatrix,
                                          std::vector<int> sequence) {
  std::vector<std::vector<double>> delta(sequence.size(), std::vector<double>(transitionMatrix.size()));
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
      double p = 0;
      int m = 0;
      for (unsigned int k = 0; k < transitionMatrix.size(); ++k) {
        double q = delta[i - 1][k] * transitionMatrix[k][j];
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
  double p = 0;
  int m = 0;
  for (unsigned int i = 0; i < transitionMatrix.size(); ++i) {
    double q = delta[sequence.size() - 1][i];
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
