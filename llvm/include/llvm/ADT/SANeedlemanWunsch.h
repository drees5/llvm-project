#ifndef NEEDMAN_WUNSCH_H
#define NEEDMAN_WUNSCH_H

#include <algorithm>
#include <cstddef>
#include <functional>
#include <optional>
#include <vector>

#include "llvm/ADT/SequenceAlignment.h"

namespace {

const auto &GetDiagonal = [](const auto &Matrix, int Row, int Col) {
  return Matrix(Row - 1, Col - 1);
};
const auto &GetLeft = [](const auto &Matrix, int Row, int Col) {
  return Matrix(Row, Col - 1);
};
const auto &GetUpper = [](const auto &Matrix, int Row, int Col) {
  return Matrix(Row - 1, Col);
};

const auto &GetScoringInfo = [](const llvm::ScoringSystem &Scoring) {
  auto [Gap, Match, Mismatch, AllowMismatch] = Scoring;

  Mismatch =
      AllowMismatch ? Mismatch : std::numeric_limits<ScoreSystemType>::min();

  return std::tuple{Gap, Match, AllowMismatch, Mismatch};
};

template <typename T> class Matrix {
public:
  Matrix(size_t Rows, size_t Cols)
      : Ts{new T[Rows * Cols]}, Rows{Rows}, Cols{Cols} {};

  T &operator()(int Row, int Col) { return Ts[Row * Cols + Col]; }
  T operator()(int Row, int Col) const { return Ts[Row * Cols + Col]; }

  size_t getRows() const { return Rows; }
  size_t getCols() const { return Cols; }

private:
  T *Ts;
  size_t Rows;
  size_t Cols;
};

} // namespace

namespace llvm {

template <typename ContainerType, typename Ty, Ty, typename MatchFnTy>
class SequenceAligner;

template <typename ContainerType,
          typename Ty = typename ContainerType::value_type, Ty Blank = Ty(0),
          typename MatchFnTy = std::function<bool(Ty, Ty)>>
class NeedlemanWunschSA
    : public SequenceAligner<ContainerType, Ty, Blank, MatchFnTy> {
private:
  using BaseType = SequenceAligner<ContainerType, Ty, Blank, MatchFnTy>;

  using ScoreMatrix = Matrix<ScoreSystemType>;
  using MatchMatrix = Matrix<bool>;

  std::optional<MatchMatrix> cacheAllMatches(ContainerType &Seq1,
                                             ContainerType &Seq2) {
    if (BaseType::getMatchOperation() == nullptr) {
      return std::nullopt;
    }

    const auto Rows = std::size(Seq1), Cols = std::size(Seq2);

    auto Table = Matrix<bool>{Rows, Cols};

    for (unsigned Row = 0; Row < Rows; Row++)
      for (unsigned Col = 0; Col < Cols; Col++) {
        Table(Row, Col) = BaseType::match(Seq1[Row], Seq2[Col]);
      }

    return Table;
  }

  ScoreMatrix
  computeScoreMatrix(ContainerType &Seq1, ContainerType &Seq2,
                     const std::optional<MatchMatrix> &PossibleMatches) {

    ScoringSystem &Scoring = BaseType::getScoring();
    const auto &[Gap, Match, AllowMismatch, Mismatch] = GetScoringInfo(Scoring);

    auto Matrix = ScoreMatrix{std::size(Seq1) + 1, std::size(Seq2) + 1};

    // First element of each row
    for (size_t Row = 0; Row < Matrix.getRows(); Row++)
      Matrix(Row, 0) = Row * Gap;
    // First row
    for (size_t Col = 0; Col < Matrix.getCols(); Col++)
      Matrix(0, Col) = Col * Gap;

    const auto &IndexDiagonalMatches = [&PossibleMatches, &Seq1,
                                        &Seq2](int Row, int Column) {
      if (PossibleMatches) {
        return GetDiagonal(*PossibleMatches, Row, Column);
      }

      return Seq1[Row - 1] == Seq2[Column - 1];
    };

    const auto &GetDiagonalScore = [AllowMismatch, Match,
                                    Mismatch](auto DiagonalMatches,
                                              auto DiagonalScore) {
      if (AllowMismatch) {
        ScoreSystemType Similarity = DiagonalMatches ? Match : Mismatch;

        return DiagonalScore + Similarity;
      }

      return DiagonalMatches ? DiagonalScore : Mismatch;
    };

    for (unsigned Row = 1; Row < Matrix.getRows(); Row++) {
      for (unsigned Col = 1; Col < Matrix.getCols(); Col++) {
        ScoreSystemType Diagonal = GetDiagonalScore(
            IndexDiagonalMatches(Row, Col), GetDiagonal(Matrix, Row, Col));
        ScoreSystemType Left = GetLeft(Matrix, Row, Col) + Gap;
        ScoreSystemType Upper = GetUpper(Matrix, Row, Col) + Gap;

        ScoreSystemType Score = std::max({Diagonal, Upper, Left});

        Matrix(Row, Col) = Score;
      }
    }

    return Matrix;
  }

  AlignedSequence<Ty, Blank>
  buildResult(ContainerType &Seq1, ContainerType &Seq2,
              const std::optional<MatchMatrix> &PossibleMatches,
              const ScoreMatrix &Scores) {
    AlignedSequence<Ty, Blank> Result;
    auto &Data = Result.Data;

    ScoringSystem &Scoring = BaseType::getScoring();
    const auto [Gap, Match, AllowMismatch, Mismatch] = GetScoringInfo(Scoring);

    int Row = Scores.getRows() - 1, Column = Scores.getCols() - 1;

    const auto &IsDiagonal = [](auto Row, auto Column) -> bool {
      return Row > 0 && Column > 0;
    };
    const auto &IsUp = [&Scores, Gap](auto Row, auto Column) -> bool {
      return Row > 0 &&
             Scores(Row, Column) == (GetUpper(Scores, Row, Column) + Gap);
    };
    const auto &IsLeft = [&Scores, Gap](auto Row, auto Column) -> bool {
      return Column > 0 &&
             Scores(Row, Column) == (GetLeft(Scores, Row, Column) + Gap);
    };

    while (Row > 0 || Column > 0) {
      if (IsDiagonal(Row, Column)) {
        bool IsValidMatch = PossibleMatches
                                ? GetDiagonal(*PossibleMatches, Row, Column)
                                : Seq1[Row - 1] == Seq2[Column - 1];

        // Update data and score
        ScoreSystemType Score =
            AllowMismatch ? GetDiagonal(Scores, Row, Column) +
                                (IsValidMatch ? Match : Mismatch)
            : IsValidMatch ? (GetDiagonal(Scores, Row, Column) + Match)
                           : Mismatch;

        if (Scores(Row, Column) == Score) {
          if (IsValidMatch || AllowMismatch) {
            Data.emplace_front(Seq1[Row - 1], Seq2[Column - 1], IsValidMatch);
          } else {
            Data.emplace_front(Seq1[Row - 1], Blank, false);
            Data.emplace_front(Blank, Seq2[Column - 1], false);
          }
        }
        Row--;
        Column--;
        continue;
      }
      if (IsUp(Row, Column)) {
        Data.emplace_front(Seq1[Row - 1], Blank, false);
        Row--;
      } else if (IsLeft(Row, Column)) {
        Data.emplace_front(Blank, Seq2[Column - 1], false);
        Column--;
      }
    }

    return Result;
  }

public:
  static ScoringSystem getDefaultScoring() { return {-1, 2, -1}; }

  NeedlemanWunschSA() : BaseType(getDefaultScoring(), nullptr) {}

  NeedlemanWunschSA(ScoringSystem Scoring, MatchFnTy Match = nullptr)
      : BaseType(Scoring, Match) {}

  virtual size_t getMemoryRequirement(ContainerType &Seq1,
                                      ContainerType &Seq2) override {
    const size_t SizeSeq1 = std::size(Seq1);
    const size_t SizeSeq2 = std::size(Seq2);
    size_t MemorySize = 0;

    MemorySize += sizeof(ScoreSystemType) * (SizeSeq1 + 1) * (SizeSeq2 + 1);

    if (BaseType::getMatchOperation() != nullptr)
      MemorySize += SizeSeq1 * SizeSeq2 * sizeof(bool);

    return MemorySize;
  }

  virtual AlignedSequence<Ty, Blank>
  getAlignment(ContainerType &Seq1, ContainerType &Seq2) override {
    const auto &PossibleMatches = cacheAllMatches(Seq1, Seq2);

    const auto &Scores = computeScoreMatrix(Seq1, Seq2, PossibleMatches);

    return buildResult(Seq1, Seq2, PossibleMatches, Scores);
  }
};

} // namespace llvm
#endif
