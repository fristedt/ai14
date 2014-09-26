#ifndef _CHECKERS_PLAYER_HPP_
#define _CHECKERS_PLAYER_HPP_

#include "constants.hpp"
#include "deadline.hpp"
#include "move.hpp"
#include "gamestate.hpp"
#include <vector>
#include <cstdlib>
#include <climits>
#include <algorithm>

namespace checkers
{

  class Player
  {
    public:
      ///perform a move
      ///\param pState the current state of the board
      ///\param pDue time before which we must have returned
      ///\return the next state the board is in after our move
      GameState play(const GameState &pState, const Deadline &pDue);
    private:

      uint8_t currentPlayer, currentOpponent;

      int eval(const GameState &state);
      GameState alphaBetaIDS(GameState pState, Deadline pDue);
      int alphabeta(const GameState &gameState, int depth, int alpha, int beta, bool maximizingPlayer);
  };


  /*namespace checkers*/ }

#endif
