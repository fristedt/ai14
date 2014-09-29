#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <cmath>

std::vector<double> getNextEmissionDist(std::vector<std::vector<double>> transitionMatrix, 
                                        std::vector<std::vector<double>> emissionMatrix, 
                                        std::vector<std::vector<double>> initialMatrix);

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

  std::vector<std::vector<double>> alpha(sequenceLength, std::vector<double>(transitionRows));
  // Initialize. Calculate alphaduder sequence[0].
  for (int i = 0; i < transitionRows; ++i) {
    alpha[0][i] = initials[0][i] * emissions[i][sequence[0]];
    // std::cout << "alphai: " << alpha[0][i] << std::endl;
  }

  for (int i = 1; i < sequenceLength; ++i) {
    for (int j = 0; j < transitionRows; ++j) {
      double sum = 0;
      for (int k = 0; k < transitionRows; ++k) {
        // std::cout << "alhpai-1k: " << alpha[i - 1][k] << std::endl;
        // std::cout << "trans    : " << transitions[k][j] << std::endl << std::endl;
        sum += alpha[i - 1][k] * transitions[k][j];
      }
      alpha[i][j] = sum * emissions[j][sequence[i]];
      // std::cout << "sum    : " << sum << std::endl;
      // std::cout << "seqi   : " << sequence[i] << std::endl;
      // std::cout << "alphaij: " << alpha[i][j] << std::endl;
    }
  }
  
  // for (int i = 0; i < sequenceLength; ++i) {
  //   for (int j = 0; j < transitionRows; ++j) {
  //     std::cout << "alpha(" << i << ", " << j << "): " << alpha[i][j];
  //     std::cout << std::endl;
  //   }
  //   std::cout << std::endl;
  // }

  double sum = 0;
  for (int i = 0; i < transitionRows; ++i) {
    sum += alpha[sequenceLength - 1][i];
  }
  std::cout << std::fixed << std::setprecision(6) << sum << std::endl;

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
