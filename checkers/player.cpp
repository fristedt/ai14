#include "player.hpp"

namespace checkers
{
  GameState Player::play(const GameState &pState,const Deadline &pDue)
  {
    return alphaBetaIDS(pState, pDue);
  }

  GameState Player::alphaBetaIDS(GameState pState, Deadline pDue) {
    std::vector<GameState> lNextStates;
    pState.findPossibleMoves(lNextStates);

    currentPlayer = pState.getNextPlayer();
    currentOpponent = (CELL_RED | CELL_WHITE) & ~currentPlayer;

    if (lNextStates.size() == 0) return GameState(pState, Move());

    int depth = 5;
    GameState bestState;
    int bestScore = INT_MIN;
    int score;

    while (pDue - pDue.now() > 500000) {
      for (unsigned int i = 0; i < lNextStates.size(); ++i) {
        score = alphabeta(lNextStates[i], depth, INT_MIN, INT_MAX, false);
        if (score > bestScore) {
          bestScore = score;
          bestState = lNextStates[i];
        }
      }
      ++depth;
    }

    return bestState;
  }

  int Player::eval(const GameState &state) {
    if ((currentPlayer & CELL_WHITE) && state.isWhiteWin()) return INT_MAX;
    if ((currentPlayer & CELL_RED) && state.isRedWin()) return INT_MAX;

    int totalScore = 0;
    for (int i = 1; i <= 32; ++i) {
      if (state.at(i) & currentPlayer) ++totalScore;
    }
    return totalScore;
  }

  int Player::alphabeta(const GameState &gameState, int depth, int alpha, int beta, bool maximizingPlayer) {
    if (depth == 0 || gameState.isEOG()) return eval(gameState);

    std::vector<GameState> lNextStates; 
    gameState.findPossibleMoves(lNextStates);

    if (maximizingPlayer) {
      for (unsigned int i = 0; i < lNextStates.size(); i++) {
        alpha = std::max(alpha, alphabeta(lNextStates[i], depth - 1, alpha, beta, !maximizingPlayer));
        if (alpha >= beta) return beta;
      }
      return alpha;
    } else {
      for (unsigned int i = 0; i < lNextStates.size(); i++) {
        beta = std::min(beta, alphabeta(lNextStates[i], depth - 1, alpha, beta, !maximizingPlayer));
        if (beta <= alpha) return alpha;
      }
      return beta;
    }
  }

  /*namespace checkers*/ }
