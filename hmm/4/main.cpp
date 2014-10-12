#include "hmm.h"

// Most of the stuff in this project is stolen from Nikolai Shokhirev.
// His code was written in Pascal though, so some things might differ.
// http://www.shokhirev.com/nikolai/abc/alg/hmm/hmm.html
int main(int argc, char *argv[]) {
  // Read stuff from cin. 
  string transitionString;
  getline(cin, transitionString);
  
  string emissionString;
  getline(cin, emissionString);

  string initialString;
  getline(cin, initialString);

  string sequenceString;
  getline(cin, sequenceString);

  // Transitions.
  int transitionRows, transitionColumns;
  istringstream iss(transitionString);
  iss >> transitionRows;
  iss >> transitionColumns;

  int x, y; 
  long double tmp;
  vector<vector<long double>> transitions(transitionRows, vector<long double>(transitionColumns));
  x = y = 0;
  while (iss >> tmp) {
    transitions[x][y++] = tmp;
    if (y == transitionColumns) { ++x; y = 0; }
  }

  // Emissions.
  int emissionRows, emissionColumns;
  istringstream iss1(emissionString);
  iss1 >> emissionRows;
  iss1 >> emissionColumns;

  vector<vector<long double>> emissions(emissionRows, vector<long double>(emissionColumns));
  x = y = 0;
  while (iss1 >> tmp) {
    emissions[x][y++] = tmp;
    if (y == emissionColumns) { ++x; y = 0; }
  }

  // Initials.
  int initialRows, initialColumns;
  istringstream iss2(initialString);
  iss2 >> initialRows;
  iss2 >> initialColumns;

  vector<vector<long double>> initials(initialRows, vector<long double>(initialColumns));
  x = y = 0;
  while (iss2 >> tmp) {
    initials[x][y++] = tmp;
    if (y == initialColumns) { ++x; y = 0; }
  }

  // Sequence.
  int sequenceLength;
  istringstream iss3(sequenceString);
  iss3 >> sequenceLength;
  vector<int> sequence(sequenceLength);
  x = 0;
  while (iss3 >> tmp) {
    sequence[x++] = tmp;
  }

  vector<vector<long double>> newTransitions(transitions.size(), vector<long double>(transitions[0].size()));
  vector<vector<long double>> newEmissions(emissions.size(), vector<long double>(emissions[0].size()));
  vector<vector<long double>> newInitials(initials.size(), vector<long double>(initials[0].size()));

  HMM hmm(transitions, emissions, initials);

  hmm.estimateMatrices(sequence);

  hmm.printHMM4();

  return 0;
}
