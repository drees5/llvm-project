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
#include "llvm/IR/Value.h"
#include "llvm/Support/raw_ostream.h"

namespace {

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
  MATCH_FUNCTION Match = nullptr;
  using AlignedSequenceData =
      std::list<typename AlignedSequence<TY, Blank>::Entry>;

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
      // We don't include the final element as it's always matching, and always
      // a branch instruction which cannot be moved
      PairIndexes.push_back({CurrentMergedBlockIndex, I - 1});
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

      // Shouldn't happen
      if (!OriginalIndex) {
        return Indexes;
      }

      // instructions between new-old to the right. Works under the assumption
      // This places our old instrucion at the new position, and shifts all
      // that NewIndex < OldIndex
      auto It = std::begin(Indexes.ToOriginalIndexes);
      std::advance(It, OldIndex + 1);
      auto ReverseIt = std::make_reverse_iterator(It);

      auto NewIt = std::begin(Indexes.ToOriginalIndexes);
      std::advance(NewIt, NewIndex);
      auto ReverseNewIt = std::make_reverse_iterator(NewIt);

      std::rotate(ReverseIt, ReverseIt + 1, ReverseNewIt);

      // Now, instead of trying to increment affected indices,
      // we recalc OriginalToNewIndexes completely.
      // For every merged position, if there is a corresponding original index,
      // update the mapping.
      for (size_t MergedIdx = NewIndex; MergedIdx <= OldIndex; MergedIdx++) {
        if (auto Orig = Indexes.ToOriginalIndexes[MergedIdx];
            Orig.has_value()) {
          Indexes.OriginalToNewIndexes[*Orig] = MergedIdx;
        }
      }

      return Indexes;
    };

    return {UpdateIndexes(std::move(Seq1Indexes)),
            UpdateIndexes(std::move(Seq2Indexes))};
  }

  // Returns the updated previous merged block and indexes, as we may have moved
  // things around
  static std::tuple<MergedBlock, BothPreMergeIndexes, bool>
  moveInstructionsToAboveMergedBlock(
      AlignedSequenceData &InstructionPairs, const MergedBlock &TargetBlock,
      MergedBlock PreviousBlock, BothPreMergeIndexes PreMergeIndexes,
      INSTRUCTION_DEPENDENCY_CALCULATOR &BB1InstructionDependenies,
      INSTRUCTION_DEPENDENCY_CALCULATOR &BB2InstructionDependenies) {

    const auto &IsValidMove = [&PreMergeIndexes, &BB1InstructionDependenies,
                               &BB2InstructionDependenies,
                               &PreviousBlock](int TargetIndex) {
      const auto OriginalIndex1 =
          PreMergeIndexes.first.ToOriginalIndexes[TargetIndex];
      const auto OriginalIndex2 =
          PreMergeIndexes.second.ToOriginalIndexes[TargetIndex];

      // We can only move InstructionPairs from matched blocks
      if (!(OriginalIndex1 && OriginalIndex2)) {
        return false;
      }

      const auto Instruction1ClosestDependency =
          BB1InstructionDependenies.getDependent(*OriginalIndex1);

      const auto Instruction2ClosestDependency =
          BB2InstructionDependenies.getDependent(*OriginalIndex2);

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
      return NewDependencyIndex1 < PreviousBlock.EndIndex &&
             NewDependencyIndex2 < PreviousBlock.EndIndex;
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

          errs() << "Moving instruction pair from " << TargetIndex << " to "
                 << PreviousBlock.EndIndex << " block info=\n";
          OldIt->dump();

          InstructionPairs.splice(NewIt, InstructionPairs, OldIt);
          return updateIndexesAfterInstructionMove(
              std::move(PreMergeIndexes), TargetIndex, PreviousBlock.EndIndex);
        };

    bool MovedAnyInstructions = false;
    for (auto I = TargetBlock.StartIndex; I < TargetBlock.EndIndex; I++) {
      if (IsValidMove(I)) {
        // Update indexes and previous block end index
        PreMergeIndexes =
            MoveToPreviousBasicBlock(I, std::move(PreMergeIndexes));
        PreviousBlock.EndIndex++;
        MovedAnyInstructions = true;
      }
    }

    return {std::move(PreviousBlock), std::move(PreMergeIndexes),
            MovedAnyInstructions};
  }

  bool reorderMergedInstructionPairs(
      AlignedSequenceData &InstructionPairs,
      INSTRUCTION_DEPENDENCY_CALCULATOR BB1InstructionDependenies,
      INSTRUCTION_DEPENDENCY_CALCULATOR BB2InstructionDependenies) {
    auto MergedBlockIndexes = findSuccessfullyMergedBlocks(InstructionPairs);

    // Early exit out if there's no merged instructions
    if (MergedBlockIndexes.empty()) {
      return false;
    }

    auto MergeIndexes = calculatePreMergeIndexes(InstructionPairs);

    bool MovedAnyInstructions = false;

    for (auto I = std::rbegin(MergedBlockIndexes),
              // We use E as std::rend as we want to stop when we get to the
              // first block, as there's no previous one
         E = std::prev(std::rend(MergedBlockIndexes));
         I != E; I++) {
      auto NextI = I;
      NextI++;

      // We need to update the prevous basic block, as we may have added
      // additional instructions to it.
      auto [UpdatedPreviousBlock, UpdatedIndexes,
            MovedAnyInstructionsInBasicBlock] =
          moveInstructionsToAboveMergedBlock(
              InstructionPairs, *I, *NextI, MergeIndexes,
              BB1InstructionDependenies, BB2InstructionDependenies);

      *NextI = std::move(UpdatedPreviousBlock);
      MergeIndexes = std::move(UpdatedIndexes);
      if (MovedAnyInstructionsInBasicBlock) {
        MovedAnyInstructions = true;
      }
    }

    return MovedAnyInstructions;
  }

public:
  MixedOperationsSequenceAligner() {}

  AlignedSequence<TY, Blank>
  getAlignment(AlignedSequence<TY, Blank> PreReorderResult,
               const CONTAINER_TYPE &Seq1, const CONTAINER_TYPE &Seq2) {

    // No need to figure out dependencies if we haven't merged anything
    if (PreReorderResult.size() == 0) {
      return PreReorderResult;
    }

    const auto &BB1InstructionDependenies =
        INSTRUCTION_DEPENDENCY_CALCULATOR{Seq1};
    const auto &BB2InstructionDependenies =
        INSTRUCTION_DEPENDENCY_CALCULATOR{Seq2};

    const auto SequenceCopy = PreReorderResult;

    if (reorderMergedInstructionPairs(PreReorderResult.Data,
                                      BB1InstructionDependenies,
                                      BB2InstructionDependenies)) {

      errs() << "Before reordering \n";
      for (const auto &Pair : SequenceCopy) {
        Pair.dump();
      }
      errs() << "After reordering \n";
      for (const auto &Pair : PreReorderResult) {
        Pair.dump();
      }
    }

    return PreReorderResult;
  }
};

#endif
} // namespace llvm
