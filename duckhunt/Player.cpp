#include "Player.hpp"
#include <cstdlib>
#include <iostream>

namespace ducks
{
  Player::Player()
  {
  }

  Action Player::shoot(const GameState &pState, const Deadline &pDue) {
    if (lastRound != pState.getRound()) {
      hmm.clear();

      for (unsigned int i = 0; i < pState.getNumBirds(); ++i) {
        HMM hmm1(this -> STATES, this -> DIRECTIONS);
        hmm.push_back(hmm1);
      }
    }

    for (unsigned int i = 0; i < pState.getNumBirds(); ++i) {
      Bird b = pState.getBird(i);
      if (b.isDead())
        continue;

      int n = b.getSeqLength();
      if (n < THRESHOLD)
        continue;

      std::vector<int> seq(n, 0);
      for (int j = 0; j < n; ++j) {
        seq[j] = b.getObservation(j);
      }

      hmm[i].estimateMatrices(seq);

      std::vector<long double> em = hmm[i].getNextEmissionDist();
      long double max = -1;
      int nextMove = -1;
      for (int j = 0; j < DIRECTIONS; ++j) {
        if (em[j] > max) {
          max = em[j];
          nextMove = j;
        }
      }

      if (max > MIN_PROBABILITY) {
        return Action(i, (EMovement)nextMove);
      }
    }

    // This line choose not to shoot
    return cDontShoot;
  }
  std::vector<ESpecies> Player::guess(const GameState &pState, const Deadline &pDue)
  {
    /*
     * Here you should write your clever algorithms to guess the species of each bird.
     * This skeleton makes no guesses, better safe than sorry!
     */

    std::vector<ESpecies> lGuesses(pState.getNumBirds(), SPECIES_UNKNOWN);
    return lGuesses;
  }

  void Player::hit(const GameState &pState, int pBird, const Deadline &pDue)
  {
    /*
     * If you hit the bird you are trying to shoot, you will be notified through this function.
     */
    std::cerr << "HIT BIRD!!!" << std::endl;
  }

  void Player::reveal(const GameState &pState, const std::vector<ESpecies> &pSpecies, const Deadline &pDue)
  {
    /*
     * If you made any guesses, you will find out the true species of those birds in this function.
     */
  }


} /*namespace ducks*/
