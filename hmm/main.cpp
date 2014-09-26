#include <iostream>
#include <sstream>
#include <vector>

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

  // Calculate distribution for all next states.
  std::vector<double> nextStateDist(transitionRows);
  for (int i = 0; i < transitionRows; ++i) {
    for (int j = 0; j < transitionColumns; ++j) {
      nextStateDist[i] += (initials[0][j] * transitions[j][i]);
    }
  }

  // Calculate distribution for all next emissions.
  std::vector<double>nextEmissionDist(emissionColumns);
  for (int i = 0; i < emissionColumns; ++i) {
    for (int j = 0; j < emissionRows; ++j) {
      nextEmissionDist[i] += nextStateDist[j] * emissions[j][i];
    }
  }

  // Print that shit.
  std::cout << 1 << " " << emissionColumns;
  for (int i = 0; i < emissionColumns; ++i) {
    std::cout << " " << nextEmissionDist[i];
  }

  return 0;
}
