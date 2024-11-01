
#include "llvm/ADT/SequenceAlignment.h"
#include <functional>
namespace llvm {

template <typename ContainerType,
          typename Ty = typename ContainerType::value_type, Ty Blank = Ty(0),
          typename MatchFnTy = std::function<bool(Ty, Ty)>>
class NeedlemanWunschSA
    : public SequenceAligner<ContainerType, Ty, Blank, MatchFnTy> {
private:
  ScoreSystemType *Matrix;
  size_t MatrixRows;
  size_t MatrixCols;
  bool *Matches;
  size_t MatchesRows;
  size_t MatchesCols;

  const static unsigned END = 0;
  const static unsigned DIAGONAL = 1;
  const static unsigned UP = 2;
  const static unsigned LEFT = 3;

  size_t MaxRow;
  size_t MaxCol;

  using BaseType = SequenceAligner<ContainerType, Ty, Blank, MatchFnTy>;

  void cacheAllMatches(ContainerType &Seq1, ContainerType &Seq2) {
    if (BaseType::getMatchOperation() == nullptr) {
      Matches = nullptr;
      return;
    }
    const size_t SizeSeq1 = Seq1.size();
    const size_t SizeSeq2 = Seq2.size();

    MatchesRows = SizeSeq1;
    MatchesCols = SizeSeq2;
    Matches = new bool[SizeSeq1 * SizeSeq2];
    for (unsigned I = 0; I < SizeSeq1; I++)
      for (unsigned J = 0; J < SizeSeq2; J++)
        Matches[I * SizeSeq2 + J] = BaseType::match(Seq1[I], Seq2[J]);
  }

  void computeScoreMatrix(ContainerType &Seq1, ContainerType &Seq2) {
    const size_t SizeSeq1 = Seq1.size();
    const size_t SizeSeq2 = Seq2.size();

    const size_t NumRows = SizeSeq1 + 1;
    const size_t NumCols = SizeSeq2 + 1;
    Matrix = new ScoreSystemType[NumRows * NumCols];
    MatrixRows = NumRows;
    MatrixCols = NumCols;

    ScoringSystem &Scoring = BaseType::getScoring();
    const ScoreSystemType Gap = Scoring.getGapPenalty();
    const ScoreSystemType Match = Scoring.getMatchProfit();
    const bool AllowMismatch = Scoring.getAllowMismatch();
    const ScoreSystemType Mismatch =
        AllowMismatch ? Scoring.getMismatchPenalty()
                      : std::numeric_limits<ScoreSystemType>::min();

    for (unsigned I = 0; I < NumRows; I++)
      Matrix[I * NumCols + 0] = I * Gap;
    for (unsigned J = 0; J < NumCols; J++)
      Matrix[0 * NumCols + J] = J * Gap;

    ScoreSystemType MaxScore = std::numeric_limits<ScoreSystemType>::min();
    if (Matches) {
      if (AllowMismatch) {
        for (unsigned I = 1; I < NumRows; I++) {
          for (unsigned J = 1; J < NumCols; J++) {
            ScoreSystemType Similarity =
                Matches[(I - 1) * MatchesCols + J - 1] ? Match : Mismatch;
            ScoreSystemType Diagonal =
                Matrix[(I - 1) * NumCols + J - 1] + Similarity;
            ScoreSystemType Upper = Matrix[(I - 1) * NumCols + J] + Gap;
            ScoreSystemType Left = Matrix[I * NumCols + J - 1] + Gap;
            ScoreSystemType Score = std::max(std::max(Diagonal, Upper), Left);
            Matrix[I * NumCols + J] = Score;
            if (Score >= MaxScore) {
              MaxScore = Score;
              MaxRow = I;
              MaxCol = J;
            }
          }
        }
      } else {
        for (unsigned I = 1; I < NumRows; I++) {
          for (unsigned J = 1; J < NumCols; J++) {
            ScoreSystemType Diagonal =
                Matches[(I - 1) * MatchesCols + J - 1]
                    ? (Matrix[(I - 1) * NumCols + J - 1] + Match)
                    : Mismatch;
            ScoreSystemType Upper = Matrix[(I - 1) * NumCols + J] + Gap;
            ScoreSystemType Left = Matrix[I * NumCols + J - 1] + Gap;
            ScoreSystemType Score = std::max(std::max(Diagonal, Upper), Left);
            Matrix[I * NumCols + J] = Score;
            if (Score >= MaxScore) {
              MaxScore = Score;
              MaxRow = I;
              MaxCol = J;
            }
          }
        }
      }
    } else {
      if (AllowMismatch) {
        for (unsigned I = 1; I < NumRows; I++) {
          for (unsigned J = 1; J < NumCols; J++) {
            ScoreSystemType Similarity =
                (Seq1[I - 1] == Seq2[J - 1]) ? Match : Mismatch;
            ScoreSystemType Diagonal =
                Matrix[(I - 1) * NumCols + J - 1] + Similarity;
            ScoreSystemType Upper = Matrix[(I - 1) * NumCols + J] + Gap;
            ScoreSystemType Left = Matrix[I * NumCols + J - 1] + Gap;
            ScoreSystemType Score = std::max(std::max(Diagonal, Upper), Left);
            Matrix[I * NumCols + J] = Score;
            if (Score >= MaxScore) {
              MaxScore = Score;
              MaxRow = I;
              MaxCol = J;
            }
          }
        }
      } else {
        for (unsigned I = 1; I < NumRows; I++) {
          for (unsigned J = 1; J < NumCols; J++) {
            ScoreSystemType Diagonal =
                (Seq1[I - 1] == Seq2[J - 1])
                    ? (Matrix[(I - 1) * NumCols + J - 1] + Match)
                    : Mismatch;
            ScoreSystemType Upper = Matrix[(I - 1) * NumCols + J] + Gap;
            ScoreSystemType Left = Matrix[I * NumCols + J - 1] + Gap;
            ScoreSystemType Score = std::max(std::max(Diagonal, Upper), Left);
            Matrix[I * NumCols + J] = Score;
            if (Score >= MaxScore) {
              MaxScore = Score;
              MaxRow = I;
              MaxCol = J;
            }
          }
        }
      }
    }
  }

  void buildResult(ContainerType &Seq1, ContainerType &Seq2,
                   AlignedSequence<Ty, Blank> &Result) {
    auto &Data = Result.Data;

    ScoringSystem &Scoring = BaseType::getScoring();
    const ScoreSystemType Gap = Scoring.getGapPenalty();
    const ScoreSystemType Match = Scoring.getMatchProfit();
    const bool AllowMismatch = Scoring.getAllowMismatch();
    const ScoreSystemType Mismatch =
        AllowMismatch ? Scoring.getMismatchPenalty()
                      : std::numeric_limits<ScoreSystemType>::min();

    int i = MatrixRows - 1, J = MatrixCols - 1;

    size_t LongestMatch = 0;
    size_t CurrentMatch = 0;

    while (i > 0 || J > 0) {
      if (i > 0 && J > 0) {
        // Diagonal

        bool IsValidMatch = false;

        ScoreSystemType Score = std::numeric_limits<ScoreSystemType>::min();
        if (Matches) {
          IsValidMatch = Matches[(i - 1) * MatchesCols + J - 1];
        } else {
          IsValidMatch = (Seq1[i - 1] == Seq2[J - 1]);
        }

        if (!IsValidMatch) {
          if (CurrentMatch > LongestMatch)
            LongestMatch = CurrentMatch;
          CurrentMatch = 0;
        } else
          CurrentMatch += 1;

        if (AllowMismatch) {
          Score = Matrix[(i - 1) * MatrixCols + J - 1] +
                  (IsValidMatch ? Match : Mismatch);
        } else {
          Score = IsValidMatch ? (Matrix[(i - 1) * MatrixCols + J - 1] + Match)
                               : Mismatch;
        }

        if (Matrix[i * MatrixCols + J] == Score) {
          if (IsValidMatch || AllowMismatch) {
            Data.push_front(typename BaseType::EntryType(
                Seq1[i - 1], Seq2[J - 1], IsValidMatch));
          } else {
            Data.push_front(
                typename BaseType::EntryType(Seq1[i - 1], Blank, false));
            Data.push_front(
                typename BaseType::EntryType(Blank, Seq2[J - 1], false));
          }

          i--;
          J--;
          continue;
        }
      }
      if (i > 0 && Matrix[i * MatrixCols + J] ==
                       (Matrix[(i - 1) * MatrixCols + J] + Gap)) {
        // Up
        Data.push_front(
            typename BaseType::EntryType(Seq1[i - 1], Blank, false));
        i--;
      } else if (J > 0 && Matrix[i * MatrixCols + J] ==
                              (Matrix[i * MatrixCols + (J - 1)] + Gap)) {
        // Left
        Data.push_front(
            typename BaseType::EntryType(Blank, Seq2[J - 1], false));
        J--;
      }
    }

    if (CurrentMatch > LongestMatch)
      LongestMatch = CurrentMatch;
  }

  void clearAll() {
    if (Matrix)
      delete[] Matrix;
    if (Matches)
      delete[] Matches;
    Matrix = nullptr;
    Matches = nullptr;
  }

public:
  static ScoringSystem getDefaultScoring() { return ScoringSystem(-1, 2, -1); }

  NeedlemanWunschSA()
      : BaseType(getDefaultScoring(), nullptr), Matrix(nullptr),
        Matches(nullptr) {}

  NeedlemanWunschSA(ScoringSystem Scoring, MatchFnTy Match = nullptr)
      : BaseType(Scoring, Match), Matrix(nullptr), Matches(nullptr) {}

  ~NeedlemanWunschSA() { clearAll(); }

  virtual size_t getMemoryRequirement(ContainerType &Seq1,
                                      ContainerType &Seq2) override {
    const size_t SizeSeq1 = Seq1.size();
    const size_t SizeSeq2 = Seq2.size();
    size_t MemorySize = 0;

    MemorySize += sizeof(ScoreSystemType) * (SizeSeq1 + 1) * (SizeSeq2 + 1);

    if (BaseType::getMatchOperation() != nullptr)
      MemorySize += SizeSeq1 * SizeSeq2 * sizeof(bool);

    return MemorySize;
  }

  virtual AlignedSequence<Ty, Blank>
  getAlignment(ContainerType &Seq1, ContainerType &Seq2) override {
    AlignedSequence<Ty, Blank> Result;
    cacheAllMatches(Seq1, Seq2);
    computeScoreMatrix(Seq1, Seq2);
    buildResult(Seq1, Seq2, Result);
    clearAll();
    return Result;
  }
};
} // namespace llvm
