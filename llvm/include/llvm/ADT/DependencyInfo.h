#include "llvm/ADT/DenseMap.h"
#include "llvm/ADT/SmallBitVector.h"
#include "llvm/IR/InstrTypes.h"
#include "llvm/IR/Instruction.h"
#include "llvm/Support/raw_ostream.h"

#define ENABLE_DEBUG

enum class Relation {
  NONE,
  ANCESTOR,
  DESCENTANT,
};

template <typename ContainerType> class DependencyInfo {
  // Using small bit vectors might be inefficient

  llvm::SmallVector<llvm::SmallBitVector> Dep;
  llvm::SmallVector<llvm::SmallBitVector> DataDep;
  size_t Size;

public:
  using value_type = ContainerType;

  DependencyInfo(const ContainerType &Seq) : Size{Seq.size()} {

#ifdef ENABLE_DEBUG
    llvm::errs() << "Sequence Size: " << Size << "\n";
#endif

    // Instruction pointer to index in the Dep arrays mapping
    const llvm::SmallDenseMap<llvm::Instruction *, size_t> IMap = [&Seq,
                                                                   this]() {
      // We form triangular 2D arrays where row > column
      // Columns represent producing instructions,
      // Rows represent consuming instructions
      // A set bit represents a dependency between the two instructions
      llvm::SmallDenseMap<llvm::Instruction *, size_t> IMap;
      for (size_t Idx = 0; Idx < Size; ++Idx) {
        // emplace back an Idx-sized bit vector
        Dep.emplace_back(Idx);
        // Same
        DataDep.emplace_back(Idx);

        // Maintain a mapping from instruction pointer to index
        if (llvm::Instruction *I = llvm::dyn_cast<llvm::Instruction>(Seq[Idx]);
            I) {
          IMap[I] = Idx;
        }
      }
      return IMap;
    }();

    const auto InstructionReadsOrWrites = [](llvm::Instruction *I) -> bool {
      if (!I)
        return false;
      if (I->mayReadFromMemory() || I->mayWriteToMemory())
        return true;
      return false;
    };

    const auto InstructionMayLeaveBasicBlock =
        [](llvm::Instruction *I) -> bool {
      return I->mayThrow() || !I->willReturn() || llvm::isa<llvm::CallBase>(I);
    };

    const auto AddInstructionDependencies =
        [&IMap, this](const llvm::Instruction *I, const llvm::Instruction *UI,
                      size_t ProducingIdx) {
          // Find index of this instruction
          if (const auto It = IMap.find(UI); It != IMap.end()) {
            const size_t UsingJdx = It->second;

            // skip if it turns out the user it's a self-use
            if (ProducingIdx == UsingJdx)
              return;

            // The producer should be before the user
            assert(ProducingIdx < UsingJdx);
            // UsingJdx depends on ProducingIdx
            Dep[UsingJdx].set(ProducingIdx);
            DataDep[UsingJdx].set(ProducingIdx);
          }
        };

    // Optimisation: count serializing instructions and basically bail
    // out if we have too many of them Rationale: we will have few
    // opportunities for moving instructions anyway, while the cost of
    // checking the dependency info will be higher
    size_t SerializingCount = 0;

    // Iterate over all instructions in the basic block
    for (size_t ProducingIdx = 0; ProducingIdx < Size - 1; ++ProducingIdx) {
      llvm::Instruction *I =
          llvm::dyn_cast<llvm::Instruction>(Seq[ProducingIdx]);
      if (I == nullptr)
        continue;

      for (const llvm::User *U : I->users()) {
        const llvm::Instruction *UI = llvm::dyn_cast<llvm::Instruction>(U);

        // skip if this user doesn't exist or it's in a different basic block
        if (!UI || UI->getParent() != I->getParent())
          return;

        AddInstructionDependencies(I, UI, ProducingIdx);
      }

      if (InstructionMayLeaveBasicBlock) {
        // Add dependencies between this and all other instructions,
        setAll(ProducingIdx);
        SerializingCount++;
        continue;
      }

      if (I->mayWriteToMemory()) {
        // Add dependencies between this and all other memory instructions,
        setAllIf(Seq, ProducingIdx, InstructionReadsOrWrites);
      }
    }

    // If there are too many serializing instructions, then connecting all
    // indirect dependencies will take ages, while the effect will be that
    // all instructions will be dependent on almost all previous instructions.
    // Let's spare us the agony and just mark everything as dependent on
    // everything else.
    // TODO: There are more elegant ways of fixing this, e.g. treating each
    // block between two serializing instructions as a separate block for
    // alignment purposes, but this would require some rewriting
    bool TooManySerializingInstructions =
        (((Size > 64) && (SerializingCount > Size / 3)) ||
         ((Size > 512) && (SerializingCount > Size / 4)) ||
         ((Size > 1024) && (SerializingCount > 256)));

    if (TooManySerializingInstructions) {
      for (size_t Idx = 1; Idx < Size; ++Idx) {
        Dep[Idx].set();
        Dep[Idx].reset(0);
      }
      return;
    }

    // Column 0 represents the BB entry (not a real instruction)
    // It's treated separately
    for (size_t Idx = 1; Idx < Size; ++Idx)
      Dep[Idx].reset(0);

    const size_t TerminatorIndex = Size - 1;

    // Connect all ancestors and descentants
    for (size_t Idx = 2; Idx < TerminatorIndex; ++Idx)
      for (const size_t Jdx : Dep[Idx].set_bits())
        for (size_t Kdx = Idx + 1; Kdx < Size; ++Kdx)
          if (Dep[Kdx][Idx])
            Dep[Kdx].set(Jdx);

    // All instructions have a control flow dependency on the terminator
    Dep[TerminatorIndex].set();
    Dep[TerminatorIndex].reset(0);
  }

  void setAllIf(const ContainerType &Seq, size_t Idx,
                std::function<bool(llvm::Instruction *)> Predicate) {
    for (size_t Jdx = 1; Jdx < Idx; ++Jdx) {
      llvm::Instruction *ID = llvm::dyn_cast<llvm::Instruction>(Seq[Jdx]);
      if (Predicate(ID))
        Dep[Idx].set(Jdx);
    }
    for (size_t Jdx = Idx + 1; Jdx < Size; ++Jdx) {
      llvm::Instruction *ID = llvm::dyn_cast<llvm::Instruction>(Seq[Jdx]);
      if (Predicate(ID))
        Dep[Jdx].set(Idx);
    }
  }

  void setAll(size_t Idx) {
    Dep[Idx].set();
    for (size_t Jdx = Idx + 1; Jdx < Size; ++Jdx) {
      Dep[Jdx].set(Idx);
    }
  }

  bool isReady(size_t Idx) { return Dep[Idx].none(); }

  bool hasDependents(size_t Idx) {
    for (size_t Idx2 = Idx + 1; Idx2 < Dep.size(); ++Idx2)
      if (Dep[Idx2][Idx])
        return true;
    return false;
  }

  // Remove dependency between producer Idx and user Jdx
  // return true if Jdx has no remaining dependencies on previous instructions
  bool removeDependency(const size_t Idx, const size_t Jdx) {
#ifdef ENABLE_DEBUG
    // llvm::errs() << "BEFOR: " << (Dep[Jdx][Idx] ? "Y" : "N") << "\t";
    // printStateRow(Idx);
#endif
    if (!Dep[Jdx][Idx])
      return false;

    Dep[Jdx].reset(Idx);
    return isReady(Jdx);
  }

  // Remove all dependencies between producer Idx and later instructions
  // Return the instructions that have now no dependencies or earlier
  // instructions Meaningful if you are scheduling instructions in program order
  llvm::SmallVector<size_t> removeDependencies(const size_t Idx) {
    llvm::SmallVector<size_t> Result;
    llvm::SmallBitVector Cpy = Dep[Idx];

#ifdef ENABLE_GA_DEBUG
    llvm::errs() << "RemoveDeps\n";
    for (size_t Jdx = Idx + 1; Jdx < Size; ++Jdx) {
      llvm::errs() << (Dep[Jdx][Idx] ? "Y " : "N ");
    }
    printStateRow(Idx);
#endif

    Dep[Idx].clear();
    for (size_t Jdx = Idx + 1; Jdx < Size; ++Jdx) {
      if (Dep[Jdx][Idx]) {
        Dep[Jdx][Idx] = false;
        if (isReady(Jdx)) {
          Result.push_back(Jdx);
        }
      }
    }

#ifdef ENABLE_GA_DEBUG
    for (size_t i : Result)
      llvm::errs() << i << " ";
    llvm::errs() << "\n";
#endif

    return Result;
  }
  // Remove all dependencies between consumer Idx and previous instructions
  // Return the instructions that have now no dependents
  // Meaningful if you are scheduling instructions in reverse order
  llvm::SmallVector<size_t> removeDependenciesReverse(const size_t Idx) {
    llvm::SmallVector<size_t> Result;
    llvm::SmallBitVector Cpy = Dep[Idx];

#ifdef ENABLE_GA_DEBUG
    llvm::errs() << "RemoveDeps\n";
    printStateRow(Idx);
#endif

    Dep[Idx].clear();
    for (size_t Jdx : Cpy.set_bits())
      if (!hasDependents(Jdx))
        Result.push_back(Jdx);

#ifdef ENABLE_GA_DEBUG
    for (size_t i : Result)
      llvm::errs() << i << " ";
    llvm::errs() << "\n";
#endif

    return Result;
  }

  // Get all the data dependencies (forward and backward) of Idx
  llvm::SmallBitVector getConnections(size_t Idx) const {
    llvm::SmallBitVector Conns(Size);
    for (size_t OtherIdx = 0; OtherIdx < Idx; ++OtherIdx) {
      Conns[OtherIdx] = DataDep[Idx][OtherIdx];
    }
    for (size_t OtherIdx = Idx + 1; OtherIdx < Size; ++OtherIdx) {
      Conns[OtherIdx] = DataDep[OtherIdx][Idx];
    }
    return Conns;
  }

  // Returns the closest instruction that Idx depends on if any
  std::optional<size_t> getDependent(size_t Idx) const {
    if (Idx == 0) {
      return std::nullopt;
    }

    for (size_t I = Idx - 1; I--;) {
      if (getRelation(I, Idx) == Relation::DESCENTANT) {
        return I;
      }
    }
    return std::nullopt;
  }

  Relation getRelation(size_t Idx, size_t Jdx) const {
    if (Idx < Jdx && Dep[Jdx][Idx])
      return Relation::ANCESTOR;
    if (Idx > Jdx && Dep[Idx][Jdx])
      return Relation::DESCENTANT;
    return Relation::NONE;
  }

  void printStateRow(const size_t Idx) {
    for (size_t Jdx = 0; Jdx < Idx; ++Jdx) {
      llvm::errs() << (Dep[Idx][Jdx] ? "Y " : "N ");
    }
    llvm::errs() << "\n";
  }

  void printState() {
    for (size_t Idx = 0; Idx < Size; ++Idx) {
      printStateRow(Idx);
    }
  }

  size_t size() { return Size; }
};
