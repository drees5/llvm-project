#ifndef MIXED_ORDER_OPERATIONS_SEQUENCE_ALIGNER
#define MIXED_ORDER_OPERATIONS_SEQUENCE_ALIGNER

#include <algorithm>
#include <cstddef>
#include <functional>
#include <iterator>
#include <optional>
#include <unordered_map>
#include <utility>
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

struct MergedBlock {
  size_t StartIndex;
  size_t EndIndex;
};

enum class Sequence { Type1, Type2 };

template <Sequence T> struct PreMergeIndexes {
  std::vector<std::optional<int>> ToOriginalIndexes;
  std::vector<int> OriginalToNewIndexes;
};

} // namespace

namespace llvm {

template <typename INSTRUCTION_DEPENDENCY_CALCULATOR,
          typename CONTAINER_TYPE =
              typename INSTRUCTION_DEPENDENCY_CALCULATOR::value_type,
          typename TY = typename CONTAINER_TYPE::value_type, TY Blank = TY(0),
          typename MATCH_FUNCTION = std::function<bool(TY, TY)>>
class MixedOperationsSequenceAligner {
private:
  ScoringSystem Scoring = getDefaultScoring();
  MATCH_FUNCTION Match = nullptr;

  using ScoreMatrix = Matrix<ScoreSystemType>;
  using MatchMatrix = Matrix<bool>;
  using AlignedSequenceData =
      std::list<typename AlignedSequence<TY, Blank>::Entry>;

  std::optional<MatchMatrix> cacheAllMatches(const CONTAINER_TYPE &Seq1,
                                             const CONTAINER_TYPE &Seq2) {
    if (!Match) {
      return std::nullopt;
    }

    const auto Rows = std::size(Seq1), Cols = std::size(Seq2);

    auto Table = MatchMatrix{Rows, Cols};

    for (unsigned Row = 0; Row < Rows; Row++)
      for (unsigned Col = 0; Col < Cols; Col++) {
        bool IsMatch = Match(Seq1[Row], Seq2[Col]);
        Table(Row, Col) = IsMatch;
      }

    return Table;
  }

  ScoreMatrix
  computeScoreMatrix(const CONTAINER_TYPE &Seq1, const CONTAINER_TYPE &Seq2,
                     const std::optional<MatchMatrix> &PossibleMatches) {

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

      return DiagonalMatches ? DiagonalScore + Match : Mismatch;
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

  std::vector<MergedBlock> static findSuccessfullyMergedBlocks(
      const AlignedSequenceData &InstructionPairs) {
    std::vector<MergedBlock> PairIndexes;

    bool CurrentlyInMergedBlock = false;
    size_t CurrentMergedBlockIndex = 0;

    size_t I = 0;

    for (const auto &Entry : InstructionPairs) {
      if (CurrentlyInMergedBlock) {
        // Currently in merged block
        if (Entry.match()) {
          I++;
          continue;
        }
        PairIndexes.push_back({CurrentMergedBlockIndex, I});
        CurrentlyInMergedBlock = false;
      } else if (Entry.match()) {
        // Starting a new block
        CurrentMergedBlockIndex = I;
        CurrentlyInMergedBlock = true;
      }

      I++;
    }

    // This means we finished with a merged block
    if (CurrentlyInMergedBlock) {
      PairIndexes.push_back({CurrentMergedBlockIndex, I});
    }

    return PairIndexes;
  }

  using BothPreMergeIndexes = std::pair<PreMergeIndexes<Sequence::Type1>,
                                        PreMergeIndexes<Sequence::Type2>>;

  BothPreMergeIndexes static calculatePreMergeIndexes(
      const AlignedSequenceData &InstructionPairs) {
    PreMergeIndexes<Sequence::Type1> Sequence1;
    PreMergeIndexes<Sequence::Type2> Sequence2;

    size_t Seq1Index = 0, Seq2Index = 0, AlignedSeqIndex = 0;

    for (const auto &Entry : InstructionPairs) {
      if (Entry.empty()) {
      } else if (!Entry.hasBlank()) {
        Sequence1.ToOriginalIndexes.push_back(Seq1Index);
        Sequence2.ToOriginalIndexes.push_back(Seq2Index);

        Sequence1.OriginalToNewIndexes.push_back(AlignedSeqIndex);
        Sequence2.OriginalToNewIndexes.push_back(AlignedSeqIndex);

        Seq1Index++;
        Seq2Index++;
      } else if (Entry.get(0) != Blank) {
        Sequence1.ToOriginalIndexes.push_back(Seq1Index);
        Sequence2.ToOriginalIndexes.push_back(std::nullopt);

        Sequence1.OriginalToNewIndexes.push_back(AlignedSeqIndex);

        Seq1Index++;
      } else {
        Sequence2.ToOriginalIndexes.push_back(Seq2Index);
        Sequence1.ToOriginalIndexes.push_back(std::nullopt);

        Sequence2.OriginalToNewIndexes.push_back(AlignedSeqIndex);

        Seq2Index++;
      }
      AlignedSeqIndex++;
    }

    return {Sequence1, Sequence2};
  }

  static BothPreMergeIndexes
  updateIndexesAfterInstructionMove(BothPreMergeIndexes OldIndexes,
                                    size_t OldIndex, size_t NewIndex) {
    auto [Seq1Indexes, Seq2Indexes] = std::move(OldIndexes);

    const auto &UpdateIndexes = [OldIndex, NewIndex](auto Indexes) {
      auto OriginalIndex = Indexes.ToOriginalIndexes[OldIndex];
      auto NewIndexInOriginalIndex = Indexes.ToOriginalIndexes[NewIndex];

      // Shouldn't happen
      if (!OriginalIndex || !NewIndexInOriginalIndex) {
        return Indexes;
      }

      auto It = std::begin(Indexes.ToOriginalIndexes);
      std::advance(It, OldIndex + 1);
      auto ReverseIt = std::make_reverse_iterator(It);

      auto NewIt = std::begin(Indexes.ToOriginalIndexes);
      std::advance(NewIt, NewIndex);
      auto ReverseNewIt = std::make_reverse_iterator(NewIt);

      // This places our old instrucion at the new position, and shifts all
      // instructions between new-old to the right. Works under the assumption
      // that NewIndex < OldIndex
      std::rotate(ReverseIt, ReverseIt + 1, ReverseNewIt);

      // We need to increment all indexes after the new one as we've inserted a
      // new instruciton
      for (auto I = *OriginalIndex; I < *NewIndexInOriginalIndex; I++) {
        Indexes.OriginalToNewIndexes[I]++;
      }

      Indexes.OriginalToNewIndexes[*OriginalIndex] = NewIndex;

      return Indexes;
    };

    return {UpdateIndexes(std::move(Seq1Indexes)),
            UpdateIndexes(std::move(Seq2Indexes))};
  }

  // Returns the updated previous merged block and indexes, as we may have moved
  // things around
  static std::pair<MergedBlock, BothPreMergeIndexes>
  moveInstructionsToAboveMergedBlock(
      AlignedSequenceData &InstructionPairs, const MergedBlock &TargetBlock,
      MergedBlock PreviousBlock, BothPreMergeIndexes PreMergeIndexes,
      INSTRUCTION_DEPENDENCY_CALCULATOR &BB1InstructionDependenies,
      INSTRUCTION_DEPENDENCY_CALCULATOR &BB2InstructionDependenies) {

    const auto &IsValidMove = [&PreMergeIndexes, &BB1InstructionDependenies,
                               &BB2InstructionDependenies,
                               &PreviousBlock](int TargetIndex) {
      const auto OriginaIndex1 =
          PreMergeIndexes.first.ToOriginalIndexes[TargetIndex];
      const auto OriginaIndex2 =
          PreMergeIndexes.second.ToOriginalIndexes[TargetIndex];

      // We can only move InstructionPairs from matched blocks
      if (!(OriginaIndex1 && OriginaIndex2)) {
        return false;
      }

      const auto Instruction1ClosestDependency =
          BB1InstructionDependenies.getDependent(*OriginaIndex1);

      const auto Instruction2ClosestDependency =
          BB2InstructionDependenies.getDependent(*OriginaIndex2);

      // We don't mind about producers of data for minimising cross-block
      // instruction data flow, we only need to move consumers.
      if (!Instruction1ClosestDependency || !Instruction2ClosestDependency) {
        return false;
      }

      const size_t NewDependencyIndex1 =
          PreMergeIndexes.first
              .OriginalToNewIndexes[*Instruction1ClosestDependency];
      const size_t NewDependencyIndex2 =
          PreMergeIndexes.second
              .OriginalToNewIndexes[*Instruction2ClosestDependency];

      // We can move our merged instruciton iff we don't depend on anything
      // after our previous merged block.
      return NewDependencyIndex1 <= PreviousBlock.EndIndex &&
             NewDependencyIndex2 <= PreviousBlock.EndIndex;
    };

    // We don't care where in the previous basic block we've moved to, only that
    // it's now in the previous one. This is because the later compilation steps
    // can re-order the instructions within one basic block.
    const auto &MoveToPreviousBasicBlock =
        [&PreviousBlock, &InstructionPairs](
            int TargetIndex, BothPreMergeIndexes PreMergeIndexes) {
          auto NewIt = std::begin(InstructionPairs);
          std::advance(NewIt, PreviousBlock.EndIndex);

          auto OldIt = std::begin(InstructionPairs);
          std::advance(OldIt, TargetIndex);

          InstructionPairs.splice(NewIt, InstructionPairs, OldIt);
          return updateIndexesAfterInstructionMove(
              std::move(PreMergeIndexes), TargetIndex, PreviousBlock.EndIndex);
        };

    for (auto I = TargetBlock.StartIndex; I < TargetBlock.EndIndex; I++) {
      if (IsValidMove(I)) {
        // Update indexes and previous block end index
        PreMergeIndexes =
            MoveToPreviousBasicBlock(I, std::move(PreMergeIndexes));
        PreviousBlock.EndIndex++;
      }
    }

    return {std::move(PreviousBlock), std::move(PreMergeIndexes)};
  }

  void reorderMergedInstructionPairs(
      AlignedSequenceData &InstructionPairs,
      INSTRUCTION_DEPENDENCY_CALCULATOR BB1InstructionDependenies,
      INSTRUCTION_DEPENDENCY_CALCULATOR BB2InstructionDependenies) {
    auto MergedBlockIndexes = findSuccessfullyMergedBlocks(InstructionPairs);

    auto MergeIndexes = calculatePreMergeIndexes(InstructionPairs);

    for (auto I = std::rbegin(MergedBlockIndexes),
              // We use E as std::rend as we want to stop when we get to the
              // first block, as there's no previous one
         E = std::rend(MergedBlockIndexes) - 1;
         I != E; I++) {
      auto NextI = I;
      std::advance(NextI, 1);

      // We need to update the prevous basic block, as we may have added
      // additional instructions to it.
      auto [UpdatedPreviousBlock, UpdatedIndexes] =
          moveInstructionsToAboveMergedBlock(
              InstructionPairs, *I, *NextI, MergeIndexes,
              BB1InstructionDependenies, BB2InstructionDependenies);

      *NextI = std::move(UpdatedPreviousBlock);
      MergeIndexes = std::move(UpdatedIndexes);
    }
  }

  AlignedSequence<TY, Blank>
  buildResult(const CONTAINER_TYPE &Seq1, const CONTAINER_TYPE &Seq2,
              const std::optional<MatchMatrix> &PossibleMatches,
              const ScoreMatrix &Scores) {
    AlignedSequence<TY, Blank> Result{};
    auto &Data = Result.Data;

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
          Row--;
          Column--;
          continue;
        }
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

  MixedOperationsSequenceAligner() {}

  MixedOperationsSequenceAligner(ScoringSystem Scoring,
                                 MATCH_FUNCTION Match = nullptr)
      : Scoring{Scoring}, Match{Match} {}

  AlignedSequence<TY, Blank> getAlignment(const CONTAINER_TYPE &Seq1,
                                          const CONTAINER_TYPE &Seq2) {
    const auto &PossibleMatches = cacheAllMatches(Seq1, Seq2);

    const auto &Scores = computeScoreMatrix(Seq1, Seq2, PossibleMatches);

    auto PreReorderResult = buildResult(Seq1, Seq2, PossibleMatches, Scores);

    // No need to figure out dependencies if we haven't merged anything
    if (PreReorderResult.size() == 0) {
      return PreReorderResult;
    }

    const auto &BB1InstructionDependenies =
        INSTRUCTION_DEPENDENCY_CALCULATOR{Seq1};
    const auto &BB2InstructionDependenies =
        INSTRUCTION_DEPENDENCY_CALCULATOR{Seq2};

    reorderMergedInstructionPairs(PreReorderResult.Data,
                                  BB1InstructionDependenies,
                                  BB2InstructionDependenies);

    return PreReorderResult;
  }
};

#endif
} // namespace llvm
