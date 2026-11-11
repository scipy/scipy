/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file HighsSymmetry.h
 * @brief Facilities for symmetry detection
 * @author Leona Gottwald
 */

#ifndef PRESOLVE_HIGHS_SYMMETRY_H_
#define PRESOLVE_HIGHS_SYMMETRY_H_

#include <algorithm>
#include <map>
#include <vector>

#include "lp_data/HighsLp.h"
#include "util/HighsDisjointSets.h"
#include "util/HighsHash.h"
#include "util/HighsType.h"

class HighsTimer;

/// class that is responsible for assigning distinct colors for each distinct
/// double value
class HighsMatrixColoring {
  using u32 = std::uint32_t;

  std::map<double, u32> colorMap;
  double tolerance;

 public:
  // initialize with exact 0.0 and 1.0, to not have differing results due tiny
  // numerical differences on those values
  HighsMatrixColoring(double tolerance)
      : colorMap({{0.0, 1}, {1.0, 2}, {-kHighsInf, 3}, {kHighsInf, 4}}),
        tolerance(tolerance) {}

  u32 color(double value) {
    // iterator points to smallest element in map which fulfills key >= value -
    // tolerance
    auto it = colorMap.lower_bound(value - tolerance);
    // check if there is no such element, or if this element has a key value +
    // tolerance in which case we create a new color and store it with the key
    // value
    if (it == colorMap.end() || it->first > value + tolerance)
      it = colorMap.emplace_hint(it, value,
                                 static_cast<u32>(colorMap.size()) + 1);
    return it->second;
  }
};

class HighsDomain;
class HighsCliqueTable;
struct HighsSymmetries;
struct StabilizerOrbitWorkspace {
  std::vector<HighsInt> orbitPartition;
  std::vector<HighsInt> orbitSize;
  std::vector<HighsInt> linkCompressionStack;
};

struct StabilizerOrbits {
  std::vector<HighsInt> orbitCols;
  std::vector<HighsInt> orbitStarts;
  std::vector<HighsInt> stabilizedCols;
  const HighsSymmetries* symmetries;

  HighsInt orbitalFixing(HighsDomain& domain) const;

  bool isStabilized(HighsInt col) const;
};

struct HighsOrbitopeMatrix {
  enum Type {
    kFull,
    kPacking,
  };
  enum RowPackingStatus : int8_t {
    kRowUndetermined = -1,
    kRowNotPacking = 0,
    kRowPacking = 1,
    kRowPackingNegated = 2,
  };
  HighsInt rowLength;
  HighsInt numRows;
  HighsInt numSetPackingRows = 0;
  HighsHashTable<HighsInt, HighsInt> columnToRow;
  std::vector<RowPackingStatus> rowIsSetPacking;
  std::vector<HighsInt> matrix;

  HighsInt& entry(HighsInt i, HighsInt j) {
    return matrix[i + static_cast<size_t>(j) * numRows];
  }

  const HighsInt& entry(HighsInt i, HighsInt j) const {
    return matrix[i + static_cast<size_t>(j) * numRows];
  }

  HighsInt& operator()(HighsInt i, HighsInt j) { return entry(i, j); }

  const HighsInt& operator()(HighsInt i, HighsInt j) const {
    return entry(i, j);
  }

  HighsInt orbitalFixing(HighsDomain& domain) const;

  void determineOrbitopeType(HighsCliqueTable& cliquetable);

  HighsInt getBranchingColumn(const std::vector<double>& colLower,
                              const std::vector<double>& colUpper,
                              HighsInt col) const;

 private:
  HighsInt orbitalFixingForFullOrbitope(const std::vector<HighsInt>& rows,
                                        HighsDomain& domain) const;

  HighsInt orbitalFixingForPackingOrbitope(const std::vector<HighsInt>& rows,
                                           HighsDomain& domain) const;

  void detectSetPackingRows(HighsCliqueTable& cliquetable, HighsInt cliqueVal);
};

struct HighsSymmetries {
  std::vector<HighsInt> permutationColumns;
  std::vector<HighsInt> permutations;
  std::vector<HighsInt> orbitPartition;
  std::vector<HighsInt> orbitSize;
  std::vector<HighsInt> columnPosition;
  std::vector<HighsInt> linkCompressionStack;
  std::vector<HighsOrbitopeMatrix> orbitopes;
  HighsHashTable<HighsInt, HighsInt> columnToOrbitope;
  HighsInt numPerms = 0;
  HighsInt numGenerators = 0;

  void clear();
  void mergeOrbits(HighsInt col1, HighsInt col2,
                   std::vector<HighsInt>& orbitPartition,
                   std::vector<HighsInt>& orbitSize,
                   std::vector<HighsInt>& linkCompressionStack);
  void mergeOrbits(HighsInt col1, HighsInt col2) {
    mergeOrbits(col1, col2, orbitPartition, orbitSize, linkCompressionStack);
  }
  HighsInt getOrbit(HighsInt col, std::vector<HighsInt>& orbitPartition,
                    std::vector<HighsInt>& linkCompressionStack);
  HighsInt getOrbit(HighsInt col) {
    return getOrbit(col, orbitPartition, linkCompressionStack);
  }

  HighsInt propagateOrbitopes(HighsDomain& domain) const;

  HighsInt getBranchingColumn(const std::vector<double>& colLower,
                              const std::vector<double>& colUpper,
                              HighsInt col) const {
    if (columnToOrbitope.size() == 0) return col;
    const HighsInt* orbitope = columnToOrbitope.find(col);
    if (!orbitope || orbitopes[*orbitope].numSetPackingRows == 0) return col;

    return orbitopes[*orbitope].getBranchingColumn(colLower, colUpper, col);
  }

  std::shared_ptr<const StabilizerOrbits> computeStabilizerOrbits(
      const HighsDomain& localdom, StabilizerOrbitWorkspace& workspace);
};

class HighsSymmetryDetection {
  using u64 = std::uint64_t;
  using u32 = std::uint32_t;

  const HighsLp* model;
  // compressed graph storage
  std::vector<HighsInt> Gstart;
  std::vector<HighsInt> Gend;
  std::vector<std::pair<HighsInt, HighsUInt>> Gedge;

  std::vector<std::pair<HighsInt, HighsUInt>> edgeBuffer;

  std::vector<HighsInt> currentPartition;
  std::vector<HighsInt> currentPartitionLinks;
  std::vector<HighsInt> vertexToCell;
  std::vector<HighsInt> vertexPosition;
  std::vector<HighsInt> vertexGroundSet;
  std::vector<HighsInt> orbitPartition;
  std::vector<HighsInt> orbitSize;

  std::vector<HighsInt> cellCreationStack;
  std::vector<HighsBool> cellInRefinementQueue;
  std::vector<HighsInt> refinementQueue;
  std::vector<HighsInt*> distinguishCands;
  std::vector<HighsInt> automorphisms;

  std::vector<HighsInt> linkCompressionStack;

  std::vector<u32> currNodeCertificate;
  std::vector<u32> firstLeaveCertificate;
  std::vector<u32> bestLeaveCertificate;
  std::vector<HighsInt> firstLeavePartition;
  std::vector<HighsInt> bestLeavePartition;

  HighsHashTable<HighsInt, u64> vertexHash;
  HighsHashTable<std::tuple<HighsInt, HighsInt, HighsUInt>> firstLeaveGraph;
  HighsHashTable<std::tuple<HighsInt, HighsInt, HighsUInt>> bestLeaveGraph;

  HighsInt firstLeavePrefixLen;
  HighsInt bestLeavePrefixLen;
  HighsInt firstPathDepth;
  HighsInt bestPathDepth;

  HighsInt numAutomorphisms;
  HighsInt numCol;
  HighsInt numRow;
  HighsInt numVertices;
  HighsInt numActiveCols;

  // node in the search tree for finding automorphisms
  struct Node {
    HighsInt stackStart;
    HighsInt certificateEnd;
    HighsInt targetCell;
    HighsInt lastDistiguished;
  };

  std::vector<Node> nodeStack;

  HighsInt getCellStart(HighsInt pos);

  void backtrack(HighsInt backtrackStackNewEnd, HighsInt backtrackStackEnd);
  void cleanupBacktrack(HighsInt cellCreationStackPos);

  void switchToNextNode(HighsInt backtrackDepth);

  bool compareCurrentGraph(
      const HighsHashTable<std::tuple<HighsInt, HighsInt, HighsUInt>>&
          otherGraph,
      HighsInt& wrongCell) const;

  void removeFixPoints();
  void initializeGroundSet();
  HighsHashTable<std::tuple<HighsInt, HighsInt, HighsUInt>> dumpCurrentGraph();
  bool mergeOrbits(HighsInt v1, HighsInt v2);
  HighsInt getOrbit(HighsInt vertex);

  void initializeHashValues();
  void clearRefinementState();
  bool isomorphicToFirstLeave();
  bool partitionRefinement();
  bool checkStoredAutomorphism(HighsInt vertex) const;
  u64 getVertexHash(HighsInt vertex);
  HighsInt selectTargetCell() const;

  bool updateCellMembership(HighsInt vertex, HighsInt cell,
                            bool markForRefinement = true);
  bool splitCell(HighsInt cell, HighsInt splitPoint);
  void markCellForRefinement(HighsInt cell);

  bool distinguishVertex(HighsInt targetCell);
  bool determineNextToDistinguish();
  void createNode();

  HighsInt cellSize(HighsInt cell) const {
    return currentPartitionLinks[cell] - cell;
  }

  bool isFromBinaryColumn(HighsInt vertex) const;

  struct ComponentData {
    HighsDisjointSets<> components;
    std::vector<HighsInt> componentStarts;
    std::vector<HighsInt> componentSets;
    std::vector<HighsInt> componentNumOrbits;
    std::vector<HighsInt> componentNumber;
    std::vector<HighsInt> permComponentStarts;
    std::vector<HighsInt> permComponents;
    std::vector<HighsInt> firstUnfixed;
    std::vector<HighsInt> numUnfixed;

    HighsInt getComponentByIndex(HighsInt compIndex) const {
      return componentNumber[compIndex];
    }
    HighsInt numComponents() const {
      return static_cast<HighsInt>(componentStarts.size()) - 1;
    }
    HighsInt componentSize(HighsInt component) const {
      return componentStarts[component + 1] - componentStarts[component];
    }

    HighsInt getVertexComponent(HighsInt vertexPosition) {
      return components.getSet(vertexPosition);
    }

    HighsInt getPermuationComponent(HighsInt permIndex) {
      return components.getSet(firstUnfixed[permIndex]);
    }
  };

  ComponentData computeComponentData(const HighsSymmetries& symmetries);

  bool isFullOrbitope(const ComponentData& componentData, HighsInt component,
                      HighsSymmetries& symmetries);

  bool recordAutomorphism(const std::vector<HighsInt>& leavePartition,
                          HighsSymmetries& symmetries, HighsInt maxPerms,
                          HighsInt pathDepth, HighsInt& backtrackDepth);

 public:
  void loadModelAsGraph(const HighsLp& model, double epsilon);

  bool initializeDetection();

  void run(HighsSymmetries& symmetries, const HighsTimer* timer = nullptr,
           double timeLimit = kHighsInf);
};

#endif
