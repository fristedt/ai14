#include "Player.hpp"
#include <cstdlib>
#include <iostream>

namespace ducks
{
  Player::Player()
  {
  }

  Action Player::shoot(const GameState &pState, const Deadline &pDue) {
    // std::cerr << "Inside shoot." << std::endl;
    if (!initialized) {
      hmm.clear();

      for (unsigned int i = 0; i < pState.getNumBirds(); ++i) {
        HMM hmm1(this -> STATES, this -> DIRECTIONS);
        hmm.push_back(hmm1);
      }
      initialized = 1;
    }

    if (pState.getBird(0).getSeqLength() == 99) 
      initialized = 0;

    for (unsigned int i = 0; i < pState.getNumBirds(); ++i) {
      // std::cerr << "Processing bird: " << i << std::endl;
      Bird b = pState.getBird(i);
      if (b.isDead())
        continue;

      int n = b.getSeqLength();
      if (n < 76)
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

      if (max > 0.6) {
        // std::cerr << "Shooting bird " << i << std::endl;
        return Action(i, (EMovement)nextMove);
      }
    }

    // This line choose not to shoot
    return cDontShoot;
  }
  // Action Player::shoot(const GameState &pState, const Deadline &pDue)
  // {
  //   if (pState.getBird(0).getSeqLength() == 50) {
  //     initialized = 0;
  //     movementHMMs.clear();
  //   }
  //
  //   // Initialize HMMs.
  //   if (!initialized) {
  //     std::cerr << "INITIALIZING HMMS." << std::endl;
  //     for (unsigned int i = 0; i < pState.getNumBirds(); ++i) {
  //       HMM hmm(STATES, DIRECTIONS);
  //       movementHMMs.push_back(hmm);
  //     }
  //     initialized = 1; 
  //     std::cerr << "INITIALIZED." << std::endl;
  //   }
  //
  //   for (unsigned int i = 0; i < pState.getNumBirds(); ++i) {
  //     Bird bird = pState.getBird(i);
  //     if (bird.isDead()) continue;
  //
  //     if (bird.getSeqLength() > THRESHOLD) {
  //       std::vector<int> sequence;
  //       for (int j = 0; j < bird.getSeqLength(); ++j) {
  //         sequence.push_back(bird.getObservation(i));
  //       }
  //
  //       movementHMMs[i].estimateMatrices(sequence);
  //
  //       std::vector<long double> nextDist = movementHMMs[i].getNextEmissionDist();
  //
  //       long double maxVal = 0;
  //       int maxArg = 0;
  //       for (unsigned int j = 0; j < nextDist.size(); ++j) {
  //         if (nextDist[j] > maxVal) {
  //           maxVal = nextDist[j];
  //           maxArg = j;
  //         }
  //       }
  //
  //       if (maxVal >= MIN_PROBABILITY) {
  //         return Action(i, (EMovement)maxArg);
  //       }
  //     }
  //   }
  //
  //   // std::cerr << "Not shooting" << std::endl;
  //   // This line choose not to shoot
  //   return cDontShoot;
  //
  //   //This line would predict that bird 0 will move right and shoot at it
  //   // return Action(0, MOVE_RIGHT);
  // }

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
