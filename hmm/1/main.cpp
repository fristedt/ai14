#include <iostream>
#include <sstream>
#include <vector>

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

  std::vector<double> nextEmissionDist = getNextEmissionDist(transitions, emissions, initials);

  // Print that shit. HMM1.
  std::cout << 1 << " " << emissions[0].size();
  for (unsigned int i = 0; i < emissions[0].size(); ++i) {
    std::cout << " " << nextEmissionDist[i];
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
