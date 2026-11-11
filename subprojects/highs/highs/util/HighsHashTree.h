/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_UTIL_HASH_TREE_H_
#define HIGHS_UTIL_HASH_TREE_H_

#include <stdexcept>

#include "util/HighsHash.h"

using std::memcpy;
using std::memmove;

template <typename K, typename V = void>
class HighsHashTree {
  using Entry = HighsHashTableEntry<K, V>;
  using ValueType =
      typename std::remove_reference<decltype(Entry().value())>::type;

  enum Type {
    kEmpty = 0,
    kListLeaf = 1,
    kInnerLeafSizeClass1 = 2,
    kInnerLeafSizeClass2 = 3,
    kInnerLeafSizeClass3 = 4,
    kInnerLeafSizeClass4 = 5,
    kBranchNode = 6,
  };

  enum Constants {
    kBitsPerLevel = 6,
    kBranchFactor = 1 << kBitsPerLevel,
    // even though we could use up to 64 bits of the hash this would require
    // additional handling in the last levels to avoid negative shift values
    // up to 9 depth levels are Ok though as up to index 8 the
    // get_hash_chunks16() function shifts right by a non-negative amount
    kMaxDepth = 9,
    kMinLeafSize = 6,
    kLeafBurstThreshold = 54,
  };

  static uint64_t compute_hash(const K& key) {
    return HighsHashHelpers::hash(key);
  }

  static uint8_t get_hash_chunk(uint64_t hash, int pos) {
    return (hash >> (64 - kBitsPerLevel - pos * kBitsPerLevel)) &
           (kBranchFactor - 1);
  }

  static uint16_t get_hash_chunks16(uint64_t hash, int pos) {
    return (hash >> (48 - pos * kBitsPerLevel));
  }

  static uint8_t get_first_chunk16(uint16_t chunks) {
    return chunks >> (16 - kBitsPerLevel);
  }

  static void set_hash_chunk(uint64_t& hash, uint64_t chunk, int chunkPos) {
    const int shiftAmount = (60 - kBitsPerLevel - chunkPos * kBitsPerLevel);
    chunk ^= (hash >> shiftAmount) & (kBranchFactor - 1);
    hash ^= chunk << shiftAmount;
  }

  struct Occupation {
    uint64_t occupation;

    Occupation() {}
    Occupation(uint64_t occupation) : occupation(occupation) {}
    operator uint64_t() const { return occupation; }

    void set(uint8_t pos) { occupation |= uint64_t{1} << (pos); }

    void flip(uint8_t pos) { occupation ^= uint64_t{1} << (pos); }

    bool test(uint8_t pos) const { return occupation & (uint64_t{1} << pos); }

    int num_set_until(uint8_t pos) const {
      return HighsHashHelpers::popcnt(occupation >> pos);
    }

    int num_set_after(uint8_t pos) const {
      return HighsHashHelpers::popcnt(occupation << (63 - (pos)));
    }

    int num_set() const { return HighsHashHelpers::popcnt(occupation); }
  };

  static constexpr int entries_to_size_class(unsigned int numEntries) {
    return 1 + unsigned(numEntries + ((kLeafBurstThreshold - kMinLeafSize) / 3 -
                                      kMinLeafSize - 1)) /
                   ((kLeafBurstThreshold - kMinLeafSize) / 3);
  }

  template <int kSizeClass>
  struct InnerLeaf {
    static constexpr int capacity() {
      return kMinLeafSize +
             (kSizeClass - 1) * (kLeafBurstThreshold - kMinLeafSize) / 3;
    }
    // the leaf stores entries the same way as inner nodes
    // but handles collisions on the occupation flag like a
    // linear probing hash table.
    // Since the occupation flag has 64 bits and we only use
    // 15 collisions should be rare and most often we won't need
    // to do a linear scan and key comparisons at all
    Occupation occupation;
    int size;
    std::array<uint64_t, capacity() + 1> hashes;
    std::array<Entry, capacity()> entries;

    InnerLeaf() : occupation(0), size(0) { hashes[0] = 0; }

    template <int kOtherSize>
    InnerLeaf(InnerLeaf<kOtherSize>&& other) {
      assert(other.size <= capacity());
      occupation = other.occupation;
      size = other.size;
      std::copy(other.hashes.cbegin(),
                std::next(other.hashes.cbegin(), size + 1), hashes.begin());
      std::move(other.entries.begin(), std::next(other.entries.begin(), size),
                entries.begin());
    }

    int get_num_entries() const { return size; }

    std::pair<ValueType*, bool> insert_entry(uint64_t fullHash, int hashPos,
                                             Entry& entry) {
      assert(size < capacity());
      uint16_t hash = get_hash_chunks16(fullHash, hashPos);
      uint8_t hashChunk = get_first_chunk16(hash);

      int pos = occupation.num_set_until(hashChunk);

      if (occupation.test(hashChunk)) {
        // since the occupation flag is set we need to start searching from
        // pos-1 and can rely on a hash chunk with the same value existing for
        // the scan
        --pos;
        while (hashes[pos] > hash) ++pos;

        if (find_key(entry.key(), hash, pos))
          return std::make_pair(&entries[pos].value(), false);

      } else {
        occupation.set(hashChunk);

        if (pos < size)
          while (hashes[pos] > hash) ++pos;
      }

      if (pos < size) move_backward(pos, size);
      entries[pos] = std::move(entry);
      hashes[pos] = hash;
      ++size;
      hashes[size] = 0;

      return std::make_pair(&entries[pos].value(), true);
    }

    ValueType* find_entry(uint64_t fullHash, int hashPos, const K& key) {
      uint16_t hash = get_hash_chunks16(fullHash, hashPos);
      uint8_t hashChunk = get_first_chunk16(hash);
      if (!occupation.test(hashChunk)) return nullptr;

      int pos = occupation.num_set_until(hashChunk) - 1;
      while (hashes[pos] > hash) ++pos;

      if (find_key(key, hash, pos)) return &entries[pos].value();

      return nullptr;
    }

    bool erase_entry(uint64_t fullHash, int hashPos, const K& key) {
      uint16_t hash = get_hash_chunks16(fullHash, hashPos);
      uint8_t hashChunk = get_first_chunk16(hash);
      if (!occupation.test(hashChunk)) return false;

      int startPos = occupation.num_set_until(hashChunk) - 1;
      while (get_first_chunk16(hashes[startPos]) > hashChunk) ++startPos;

      int pos = startPos;
      while (hashes[pos] > hash) ++pos;

      if (!find_key(key, hash, pos)) return false;

      --size;
      if (pos < size) {
        std::move(std::next(entries.begin(), pos + 1),
                  std::next(entries.begin(), size + 1),
                  std::next(entries.begin(), pos));
        memmove(&hashes[pos], &hashes[pos + 1],
                sizeof(hashes[0]) * (size - pos));
        if (get_first_chunk16(hashes[startPos]) != hashChunk)
          occupation.flip(hashChunk);
      } else if (startPos == pos)
        occupation.flip(hashChunk);

      hashes[size] = 0;
      return true;
    }

    void rehash(int hashPos) {
      // function needs to possibly reorder elements by a different hash value
      // chances are very high we are already ordered correctly as we use 16
      // bits of the hash and one level is uses 6 bits, so the new values
      // are guaranteed to be ordered correctly by their 10 most significant
      // bits if increasing the hash position by 1 and only if the 10 bits of
      // the hash had a collision the new 6 bits might break a tie differently.
      // It is, however, important to maintain the exact ordering as otherwise
      // elements may not be found.
      occupation = 0;
      for (int i = 0; i < size; ++i) {
        hashes[i] = get_hash_chunks16(compute_hash(entries[i].key()), hashPos);
        occupation.set(get_first_chunk16(hashes[i]));
      }

      int i = 0;
      while (i < size) {
        uint8_t hashChunk = get_first_chunk16(hashes[i]);
        int pos = occupation.num_set_until(hashChunk) - 1;

        // if the position is after i the element definitely comes later, so we
        // swap it to that position and proceed without increasing i until
        // eventually an element appears that comes at position i or before
        if (pos > i) {
          std::swap(hashes[pos], hashes[i]);
          std::swap(entries[pos], entries[i]);
          continue;
        }

        // the position is before or at i, now check where the exact location
        // should be for the ordering by hash so that the invariant is that all
        // elements up to i are properly sorted. Essentially insertion sort but
        // with the modification of having a high chance to guess the correct
        // position already using the occupation flags.
        while (pos < i && hashes[pos] >= hashes[i]) ++pos;

        // if the final position is before i we need to move elements to
        // make space at that position, otherwise nothing needs to be done but
        // incrementing i increasing the sorted range by 1.
        if (pos < i) {
          uint64_t hash = hashes[i];
          auto entry = std::move(entries[i]);
          move_backward(pos, i);
          hashes[pos] = hash;
          entries[pos] = std::move(entry);
        }
        ++i;
      }
    }

    void move_backward(const int& first, const int& last) {
      // move elements backwards
      std::move_backward(std::next(entries.begin(), first),
                         std::next(entries.begin(), last),
                         std::next(entries.begin(), last + 1));
      memmove(&hashes[first + 1], &hashes[first],
              sizeof(hashes[0]) * (last - first));
    }

    bool find_key(const K& key, const uint16_t& hash, int& pos) const {
      // find key
      while (pos != size && hashes[pos] == hash) {
        if (key == entries[pos].key()) return true;
        ++pos;
      }
      return false;
    }
  };

  struct ListNode {
    ListNode* next;
    HighsHashTableEntry<K, V> entry;
    ListNode(HighsHashTableEntry<K, V>&& entry)
        : next(nullptr), entry(std::move(entry)) {}
  };
  struct ListLeaf {
    ListNode first;
    int count;

    ListLeaf(HighsHashTableEntry<K, V>&& entry)
        : first(std::move(entry)), count(1) {}
  };

  struct BranchNode;

  struct NodePtr {
    uintptr_t ptrAndType;

    NodePtr() : ptrAndType(kEmpty) {}
    NodePtr(std::nullptr_t) : ptrAndType(kEmpty) {}
    NodePtr(ListLeaf* ptr)
        : ptrAndType(reinterpret_cast<uintptr_t>(ptr) | kListLeaf) {}
    NodePtr(InnerLeaf<1>* ptr)
        : ptrAndType(reinterpret_cast<uintptr_t>(ptr) | kInnerLeafSizeClass1) {}
    NodePtr(InnerLeaf<2>* ptr)
        : ptrAndType(reinterpret_cast<uintptr_t>(ptr) | kInnerLeafSizeClass2) {}
    NodePtr(InnerLeaf<3>* ptr)
        : ptrAndType(reinterpret_cast<uintptr_t>(ptr) | kInnerLeafSizeClass3) {}
    NodePtr(InnerLeaf<4>* ptr)
        : ptrAndType(reinterpret_cast<uintptr_t>(ptr) | kInnerLeafSizeClass4) {}
    NodePtr(BranchNode* ptr)
        : ptrAndType(reinterpret_cast<uintptr_t>(ptr) | kBranchNode) {
      assert(ptr != nullptr);
    }

    Type getType() const { return Type(ptrAndType & 7u); }

    int numEntriesEstimate() const {
      switch (getType()) {
        case kEmpty:
          return 0;
        case kListLeaf:
          return 1;
        case kInnerLeafSizeClass1:
          return InnerLeaf<1>::capacity();
        case kInnerLeafSizeClass2:
          return InnerLeaf<2>::capacity();
        case kInnerLeafSizeClass3:
          return InnerLeaf<3>::capacity();
        case kInnerLeafSizeClass4:
          return InnerLeaf<4>::capacity();
        case kBranchNode:
          // something large should be returned so that the number of entries
          // is estimated above the threshold to merge when the parent checks
          // its children after deletion
          return kBranchFactor;
        default:
          throw std::logic_error("Unexpected type in hash tree");
      }
    }

    int numEntries() const {
      switch (getType()) {
        case kEmpty:
          return 0;
        case kListLeaf:
          return getListLeaf()->count;
        case kInnerLeafSizeClass1:
          return getInnerLeafSizeClass1()->size;
        case kInnerLeafSizeClass2:
          return getInnerLeafSizeClass2()->size;
        case kInnerLeafSizeClass3:
          return getInnerLeafSizeClass3()->size;
        case kInnerLeafSizeClass4:
          return getInnerLeafSizeClass4()->size;
        case kBranchNode:
          // something large should be returned so that the number of entries
          // is estimated above the threshold to merge when the parent checks
          // its children after deletion
          return kBranchFactor;
        default:
          throw std::logic_error("Unexpected type in hash tree");
      }
    }

    ListLeaf* getListLeaf() const {
      assert(getType() == kListLeaf);
      return reinterpret_cast<ListLeaf*>(ptrAndType & ~uintptr_t{7});
    }

    InnerLeaf<1>* getInnerLeafSizeClass1() const {
      assert(getType() == kInnerLeafSizeClass1);
      return reinterpret_cast<InnerLeaf<1>*>(ptrAndType & ~uintptr_t{7});
    }
    InnerLeaf<2>* getInnerLeafSizeClass2() const {
      assert(getType() == kInnerLeafSizeClass2);
      return reinterpret_cast<InnerLeaf<2>*>(ptrAndType & ~uintptr_t{7});
    }

    InnerLeaf<3>* getInnerLeafSizeClass3() const {
      assert(getType() == kInnerLeafSizeClass3);
      return reinterpret_cast<InnerLeaf<3>*>(ptrAndType & ~uintptr_t{7});
    }

    InnerLeaf<4>* getInnerLeafSizeClass4() const {
      assert(getType() == kInnerLeafSizeClass4);
      return reinterpret_cast<InnerLeaf<4>*>(ptrAndType & ~uintptr_t{7});
    }

    BranchNode* getBranchNode() const {
      assert(getType() == kBranchNode);
      return reinterpret_cast<BranchNode*>(ptrAndType & ~uintptr_t{7});
    }
  };

  struct BranchNode {
    Occupation occupation;

    NodePtr* children() { return reinterpret_cast<NodePtr*>(this + 1); }

    const NodePtr* children() const {
      return reinterpret_cast<const NodePtr*>(this + 1);
    }

    NodePtr& child(int index) { return children()[index]; }

    const NodePtr& child(int index) const { return children()[index]; }

    NodePtr* childPtr(int index) { return children() + index; }

    const NodePtr* childPtr(int index) const { return children() + index; }
  };

  static_assert(sizeof(BranchNode) % alignof(NodePtr) == 0,
                "BranchNode trailing storage must stay NodePtr aligned");

  // allocate branch nodes in multiples of 64 bytes to reduce allocator stress
  // with different sizes and reduce reallocations of nodes
  static constexpr size_t getBranchNodeSize(int numChilds) {
    return (sizeof(BranchNode) + size_t(numChilds) * sizeof(NodePtr) + 63) &
           ~size_t(63);
  };

  static BranchNode* createBranchingNode(int numChilds) {
    BranchNode* branch =
        (BranchNode*)::operator new(getBranchNodeSize(numChilds));
    branch->occupation = 0;
    return branch;
  }

  static void destroyBranchingNode(void* innerNode) {
    ::operator delete(innerNode);
  }

  static BranchNode* addChildToBranchNode(BranchNode* branch, uint8_t hashValue,
                                          int location) {
    int rightChilds = branch->occupation.num_set_after(hashValue);
    assert(rightChilds + location == branch->occupation.num_set());

    size_t newSize = getBranchNodeSize(location + rightChilds + 1);
    size_t rightSize = rightChilds * sizeof(NodePtr);

    if (newSize == getBranchNodeSize(location + rightChilds)) {
      memmove(branch->childPtr(location + 1), branch->childPtr(location),
              rightSize);

      return branch;
    }

    BranchNode* newBranch = (BranchNode*)::operator new(newSize);
    // copy the header plus the pointers left of the insertion index into the
    // new storage
    size_t leftSize = sizeof(BranchNode) + size_t(location) * sizeof(NodePtr);

    memcpy(newBranch, branch, leftSize);
    memcpy(newBranch->childPtr(location + 1), branch->childPtr(location),
           rightSize);

    destroyBranchingNode(branch);

    return newBranch;
  }

  template <int SizeClass1, int SizeClass2>
  static void mergeIntoLeaf(InnerLeaf<SizeClass1>* leaf,
                            InnerLeaf<SizeClass2>* mergeLeaf, int hashPos) {
    for (int i = 0; i < mergeLeaf->size; ++i)
      leaf->insert_entry(compute_hash(mergeLeaf->entries[i].key()), hashPos,
                         mergeLeaf->entries[i]);
  }

  template <int SizeClass>
  static void mergeIntoLeaf(InnerLeaf<SizeClass>* leaf, int hashPos,
                            NodePtr mergeNode) {
    switch (mergeNode.getType()) {
      case kListLeaf: {
        ListLeaf* mergeLeaf = mergeNode.getListLeaf();
        leaf->insert_entry(compute_hash(mergeLeaf->first.entry.key()), hashPos,
                           mergeLeaf->first.entry);
        ListNode* iter = mergeLeaf->first.next;
        while (iter != nullptr) {
          ListNode* next = iter->next;
          leaf->insert_entry(compute_hash(iter->entry.key()), hashPos,
                             iter->entry);
          delete iter;
          iter = next;
        }
        break;
      }
      case kInnerLeafSizeClass1:
        mergeIntoLeaf(leaf, mergeNode.getInnerLeafSizeClass1(), hashPos);
        delete mergeNode.getInnerLeafSizeClass1();
        break;
      case kInnerLeafSizeClass2:
        mergeIntoLeaf(leaf, mergeNode.getInnerLeafSizeClass2(), hashPos);
        delete mergeNode.getInnerLeafSizeClass2();
        break;
      case kInnerLeafSizeClass3:
        mergeIntoLeaf(leaf, mergeNode.getInnerLeafSizeClass3(), hashPos);
        delete mergeNode.getInnerLeafSizeClass3();
        break;
      case kInnerLeafSizeClass4:
        mergeIntoLeaf(leaf, mergeNode.getInnerLeafSizeClass4(), hashPos);
        delete mergeNode.getInnerLeafSizeClass4();
        break;
      default:
        break;
    }
  }

  template <int SizeClass>
  inline static bool pastChunkEnd(InnerLeaf<SizeClass>* leaf, int i, int pos) {
    return i >= leaf->size || get_first_chunk16(leaf->hashes[i]) != pos;
  }

  template <int SizeClass>
  inline static void advancePastCurrentHash(InnerLeaf<SizeClass>* leaf, int& i,
                                            int pos) {
    uint16_t commonHash = leaf->hashes[i];
    i++;
    while (!pastChunkEnd(leaf, i, pos) && leaf->hashes[i] == commonHash) ++i;
  }

  template <int SizeClass1, int SizeClass2>
  static HighsHashTableEntry<K, V>* findCommonInLeaf(
      InnerLeaf<SizeClass1>* leaf1, InnerLeaf<SizeClass2>* leaf2, int hashPos) {
    uint64_t matchMask = leaf1->occupation & leaf2->occupation;
    if (matchMask == 0) return nullptr;

    int offset1 = -1;
    int offset2 = -1;
    while (matchMask) {
      int pos = HighsHashHelpers::log2i(matchMask);
      matchMask ^= (uint64_t{1} << pos);

      int i =
          leaf1->occupation.num_set_until(static_cast<uint8_t>(pos)) + offset1;
      while (get_first_chunk16(leaf1->hashes[i]) != pos) {
        ++i;
        ++offset1;
      }

      int j =
          leaf2->occupation.num_set_until(static_cast<uint8_t>(pos)) + offset2;
      while (get_first_chunk16(leaf2->hashes[j]) != pos) {
        ++j;
        ++offset2;
      }

      while (true) {
        if (leaf1->hashes[i] > leaf2->hashes[j]) {
          ++i;
        } else if (leaf2->hashes[j] > leaf1->hashes[i]) {
          ++j;
        } else {
          int iSave = i;
          int jSave = j;
          advancePastCurrentHash(leaf1, i, pos);
          advancePastCurrentHash(leaf2, j, pos);
          for (int ii = iSave; ii < i; ++ii)
            for (int jj = jSave; jj < j; ++jj)
              if (leaf1->entries[ii].key() == leaf2->entries[jj].key())
                return &leaf1->entries[ii];
        }
        if (pastChunkEnd(leaf1, i, pos) || pastChunkEnd(leaf2, j, pos)) break;
      };
    }

    return nullptr;
  }

  template <int SizeClass>
  static HighsHashTableEntry<K, V>* findCommonInLeaf(InnerLeaf<SizeClass>* leaf,
                                                     NodePtr n2, int hashPos) {
    switch (n2.getType()) {
      case kInnerLeafSizeClass1:
        return findCommonInLeaf(leaf, n2.getInnerLeafSizeClass1(), hashPos);
      case kInnerLeafSizeClass2:
        return findCommonInLeaf(leaf, n2.getInnerLeafSizeClass2(), hashPos);
      case kInnerLeafSizeClass3:
        return findCommonInLeaf(leaf, n2.getInnerLeafSizeClass3(), hashPos);
      case kInnerLeafSizeClass4:
        return findCommonInLeaf(leaf, n2.getInnerLeafSizeClass4(), hashPos);
      case kBranchNode: {
        BranchNode* branch = n2.getBranchNode();
        uint64_t matchMask = branch->occupation & leaf->occupation;

        int offset = -1;
        while (matchMask) {
          int pos = HighsHashHelpers::log2i(matchMask);
          matchMask ^= (uint64_t{1} << pos);

          int i = leaf->occupation.num_set_until(static_cast<uint8_t>(pos)) +
                  offset;
          while (get_first_chunk16(leaf->hashes[i]) != pos) {
            ++i;
            ++offset;
          }

          int j =
              branch->occupation.num_set_until(static_cast<uint8_t>(pos)) - 1;

          do {
            if (find_recurse(branch->child(j),
                             compute_hash(leaf->entries[i].key()), hashPos + 1,
                             leaf->entries[i].key()))
              return &leaf->entries[i];
            ++i;
          } while (!pastChunkEnd(leaf, i, pos));
        }
        break;
      }
      default:
        break;
    }

    return nullptr;
  }

  static NodePtr removeChildFromBranchNode(BranchNode* branch, int location,
                                           uint64_t hash, int hashPos) {
    NodePtr newNode;
    int newNumChild = branch->occupation.num_set();

    // first check if we might be able to merge all children into one new leaf
    // based on the node numbers and assuming all of them might be in the
    // smallest size class
    if (newNumChild * InnerLeaf<1>::capacity() <= kLeafBurstThreshold) {
      // since we have a good chance of merging we now check the actual size
      // classes to see if that yields a number of entries at most the burst
      // threshold
      int childEntries = 0;
      for (int i = 0; i <= newNumChild; ++i) {
        childEntries += branch->child(i).numEntriesEstimate();
        if (childEntries > kLeafBurstThreshold) break;
      }

      if (childEntries < kLeafBurstThreshold) {
        // create a new merged inner leaf node containing all entries of
        // children first recompute the number of entries, but this time access
        // each child to get the actual number of entries needed and determine
        // this nodes size class since before we estimated the number of child
        // entries from the capacities of our child leaf node types which are
        // stored in the branch nodes pointers directly and avoid unnecessary
        // accesses of nodes that are not in cache.
        childEntries = 0;
        for (int i = 0; i <= newNumChild; ++i)
          childEntries += branch->child(i).numEntries();

        // check again if we exceed due to the extremely unlikely case
        // of having less than 5 list nodes with together more than 30 entries
        // as list nodes are only created in the last depth level
        if (childEntries < kLeafBurstThreshold) {
          switch (entries_to_size_class(childEntries)) {
            case 1: {
              InnerLeaf<1>* newLeafSize1 = new InnerLeaf<1>;
              newNode = newLeafSize1;
              for (int i = 0; i <= newNumChild; ++i)
                mergeIntoLeaf(newLeafSize1, hashPos, branch->child(i));
              break;
            }
            case 2: {
              InnerLeaf<2>* newLeafSize2 = new InnerLeaf<2>;
              newNode = newLeafSize2;
              for (int i = 0; i <= newNumChild; ++i)
                mergeIntoLeaf(newLeafSize2, hashPos, branch->child(i));
              break;
            }
            case 3: {
              InnerLeaf<3>* newLeafSize3 = new InnerLeaf<3>;
              newNode = newLeafSize3;
              for (int i = 0; i <= newNumChild; ++i)
                mergeIntoLeaf(newLeafSize3, hashPos, branch->child(i));
              break;
            }
            case 4: {
              InnerLeaf<4>* newLeafSize4 = new InnerLeaf<4>;
              newNode = newLeafSize4;
              for (int i = 0; i <= newNumChild; ++i)
                mergeIntoLeaf(newLeafSize4, hashPos, branch->child(i));
              break;
            }
            default:
              // Unexpected result from 'entries_to_size_class'
              assert(false);
          }

          destroyBranchingNode(branch);
          return newNode;
        }
      }
    }

    size_t newSize = getBranchNodeSize(newNumChild);
    size_t rightSize = (newNumChild - location) * sizeof(NodePtr);
    if (newSize == getBranchNodeSize(newNumChild + 1)) {
      // allocated size class is the same, so we do not allocate a new node
      memmove(branch->childPtr(location), branch->childPtr(location + 1),
              rightSize);
      newNode = branch;
    } else {
      // allocated size class changed, so we allocate a smaller branch node
      BranchNode* compressedBranch = (BranchNode*)::operator new(newSize);
      newNode = compressedBranch;

      size_t leftSize = sizeof(BranchNode) + size_t(location) * sizeof(NodePtr);
      memcpy(compressedBranch, branch, leftSize);
      memcpy(compressedBranch->childPtr(location),
             branch->childPtr(location + 1), rightSize);

      destroyBranchingNode(branch);
    }

    return newNode;
  }

  NodePtr root;

  template <int SizeClass>
  static std::pair<ValueType*, bool> insert_into_leaf(
      NodePtr* insertNode, InnerLeaf<SizeClass>* leaf, uint64_t hash,
      int hashPos, HighsHashTableEntry<K, V>& entry) {
    if (leaf->size == InnerLeaf<SizeClass>::capacity()) {
      auto existingEntry = leaf->find_entry(hash, hashPos, entry.key());
      if (existingEntry) return std::make_pair(existingEntry, false);

      InnerLeaf<SizeClass + 1>* newLeaf =
          new InnerLeaf<SizeClass + 1>(std::move(*leaf));
      *insertNode = newLeaf;
      delete leaf;
      return newLeaf->insert_entry(hash, hashPos, entry);
    }

    return leaf->insert_entry(hash, hashPos, entry);
  }

  static std::pair<ValueType*, bool> insert_recurse(
      NodePtr* insertNode, uint64_t hash, int hashPos,
      HighsHashTableEntry<K, V>& entry) {
    switch (insertNode->getType()) {
      case kEmpty: {
        if (hashPos == kMaxDepth) {
          ListLeaf* leaf = new ListLeaf(std::move(entry));
          *insertNode = leaf;
          return std::make_pair(&leaf->first.entry.value(), true);
        } else {
          InnerLeaf<1>* leaf = new InnerLeaf<1>;
          *insertNode = leaf;
          return leaf->insert_entry(hash, hashPos, entry);
        }
      }
      case kListLeaf: {
        ListLeaf* leaf = insertNode->getListLeaf();
        ListNode* iter = &leaf->first;
        while (true) {
          // check for existing key
          if (iter->entry.key() == entry.key())
            return std::make_pair(&iter->entry.value(), false);

          if (iter->next == nullptr) {
            // reached the end of the list and key is not duplicate, so insert
            iter->next = new ListNode(std::move(entry));
            ++leaf->count;
            return std::make_pair(&iter->next->entry.value(), true);
          }
          iter = iter->next;
        }

        break;
      }
      case kInnerLeafSizeClass1:
        return insert_into_leaf(insertNode,
                                insertNode->getInnerLeafSizeClass1(), hash,
                                hashPos, entry);
        break;
      case kInnerLeafSizeClass2:
        return insert_into_leaf(insertNode,
                                insertNode->getInnerLeafSizeClass2(), hash,
                                hashPos, entry);
        break;
      case kInnerLeafSizeClass3:
        return insert_into_leaf(insertNode,
                                insertNode->getInnerLeafSizeClass3(), hash,
                                hashPos, entry);
        break;
      case kInnerLeafSizeClass4: {
        InnerLeaf<4>* leaf = insertNode->getInnerLeafSizeClass4();
        if (leaf->size < InnerLeaf<4>::capacity())
          return leaf->insert_entry(hash, hashPos, entry);

        auto existingEntry = leaf->find_entry(hash, hashPos, entry.key());
        if (existingEntry) return std::make_pair(existingEntry, false);
        Occupation occupation = leaf->occupation;

        uint8_t hashChunk = get_hash_chunk(hash, hashPos);
        occupation.set(hashChunk);

        int branchSize = occupation.num_set();

        BranchNode* branch = createBranchingNode(branchSize);
        *insertNode = branch;
        branch->occupation = occupation;

        if (hashPos + 1 == kMaxDepth) {
          for (int i = 0; i < branchSize; ++i) branch->child(i) = nullptr;

          for (int i = 0; i < leaf->size; ++i) {
            int pos =
                occupation.num_set_until(get_first_chunk16(leaf->hashes[i])) -
                1;
            if (branch->child(pos).getType() == kEmpty)
              branch->child(pos) = new ListLeaf(std::move(leaf->entries[i]));
            else {
              ListLeaf* listLeaf = branch->child(pos).getListLeaf();
              ListNode* newNode = new ListNode(std::move(listLeaf->first));
              listLeaf->first.next = newNode;
              listLeaf->first.entry = std::move(leaf->entries[i]);
              ++listLeaf->count;
            }
          }

          delete leaf;

          ListLeaf* listLeaf;

          int pos = occupation.num_set_until(get_hash_chunk(hash, hashPos)) - 1;
          if (branch->child(pos).getType() == kEmpty) {
            listLeaf = new ListLeaf(std::move(entry));
            branch->child(pos) = listLeaf;
          } else {
            listLeaf = branch->child(pos).getListLeaf();
            ListNode* newNode = new ListNode(std::move(listLeaf->first));
            listLeaf->first.next = newNode;
            listLeaf->first.entry = std::move(entry);
            ++listLeaf->count;
          }

          return std::make_pair(&listLeaf->first.entry.value(), true);
        }

        if (branchSize > 1) {
          // maxsize in one bucket = number of items - (num buckets-1)
          // since each bucket has at least 1 item the largest one can only
          // have all remaining ones After adding the item: If it does not
          // collide
          int maxEntriesPerLeaf = 2 + leaf->size - branchSize;

          if (maxEntriesPerLeaf <= InnerLeaf<1>::capacity()) {
            // all items can go into the smallest leaf size
            for (int i = 0; i < branchSize; ++i)
              branch->child(i) = new InnerLeaf<1>;

            for (int i = 0; i < leaf->size; ++i) {
              int pos =
                  occupation.num_set_until(get_first_chunk16(leaf->hashes[i])) -
                  1;
              branch->child(pos).getInnerLeafSizeClass1()->insert_entry(
                  compute_hash(leaf->entries[i].key()), hashPos + 1,
                  leaf->entries[i]);
            }

            delete leaf;

            int pos =
                occupation.num_set_until(get_hash_chunk(hash, hashPos)) - 1;
            return branch->child(pos).getInnerLeafSizeClass1()->insert_entry(
                hash, hashPos + 1, entry);
          } else {
            // there are many collisions, determine the exact sizes first
            std::array<uint8_t, InnerLeaf<4>::capacity() + 1> sizes = {};
            sizes[occupation.num_set_until(hashChunk) - 1] += 1;
            for (int i = 0; i < leaf->size; ++i) {
              int pos =
                  occupation.num_set_until(get_first_chunk16(leaf->hashes[i])) -
                  1;
              sizes[pos] += 1;
            }

            for (int i = 0; i < branchSize; ++i) {
              switch (entries_to_size_class(sizes[i])) {
                case 1:
                  branch->child(i) = new InnerLeaf<1>;
                  break;
                case 2:
                  branch->child(i) = new InnerLeaf<2>;
                  break;
                case 3:
                  branch->child(i) = new InnerLeaf<3>;
                  break;
                case 4:
                  branch->child(i) = new InnerLeaf<4>;
                  break;
                default:
                  // Unexpected result from 'entries_to_size_class'
                  assert(false);
              }
            }

            for (int i = 0; i < leaf->size; ++i) {
              int pos =
                  occupation.num_set_until(get_first_chunk16(leaf->hashes[i])) -
                  1;

              switch (branch->child(pos).getType()) {
                case kInnerLeafSizeClass1:
                  branch->child(pos).getInnerLeafSizeClass1()->insert_entry(
                      compute_hash(leaf->entries[i].key()), hashPos + 1,
                      leaf->entries[i]);
                  break;
                case kInnerLeafSizeClass2:
                  branch->child(pos).getInnerLeafSizeClass2()->insert_entry(
                      compute_hash(leaf->entries[i].key()), hashPos + 1,
                      leaf->entries[i]);
                  break;
                case kInnerLeafSizeClass3:
                  branch->child(pos).getInnerLeafSizeClass3()->insert_entry(
                      compute_hash(leaf->entries[i].key()), hashPos + 1,
                      leaf->entries[i]);
                  break;
                case kInnerLeafSizeClass4:
                  branch->child(pos).getInnerLeafSizeClass4()->insert_entry(
                      compute_hash(leaf->entries[i].key()), hashPos + 1,
                      leaf->entries[i]);
                  break;
                default:
                  break;
              }
            }
          }

          delete leaf;

          int pos = occupation.num_set_until(hashChunk) - 1;
          insertNode = branch->childPtr(pos);
          ++hashPos;
        } else {
          // extremely unlikely that the new branch node only gets one
          // child in that case create it and defer the insertion into
          // the next depth
          branch->child(0) = leaf;
          insertNode = branch->childPtr(0);
          ++hashPos;
          leaf->rehash(hashPos);
        }

        break;
      }
      case kBranchNode: {
        BranchNode* branch = insertNode->getBranchNode();

        int location =
            branch->occupation.num_set_until(get_hash_chunk(hash, hashPos));

        if (branch->occupation.test(get_hash_chunk(hash, hashPos))) {
          --location;
        } else {
          branch = addChildToBranchNode(branch, get_hash_chunk(hash, hashPos),
                                        location);

          branch->child(location) = nullptr;
          branch->occupation.set(get_hash_chunk(hash, hashPos));
        }

        *insertNode = branch;
        insertNode = branch->childPtr(location);
        ++hashPos;
      }
    }

    return insert_recurse(insertNode, hash, hashPos, entry);
  }

  static void erase_recurse(NodePtr* erase_node, uint64_t hash, int hashPos,
                            const K& key) {
    switch (erase_node->getType()) {
      case kEmpty: {
        return;
      }
      case kListLeaf: {
        ListLeaf* leaf = erase_node->getListLeaf();

        // check for existing key
        ListNode* iter = &leaf->first;

        do {
          ListNode* next = iter->next;
          if (iter->entry.key() == key) {
            // key found, decrease count
            --leaf->count;
            if (next != nullptr) {
              // if we have a next node after replace the entry in iter by
              // moving that node into it
              *iter = std::move(*next);
              // delete memory of that node
              delete next;
            }

            break;
          }

          iter = next;
        } while (iter != nullptr);

        if (leaf->count == 0) {
          delete leaf;
          *erase_node = nullptr;
        }

        return;
      }
      case kInnerLeafSizeClass1: {
        InnerLeaf<1>* leaf = erase_node->getInnerLeafSizeClass1();
        if (leaf->erase_entry(hash, hashPos, key)) {
          if (leaf->size == 0) {
            delete leaf;
            *erase_node = nullptr;
          }
        }

        return;
      }
      case kInnerLeafSizeClass2: {
        InnerLeaf<2>* leaf = erase_node->getInnerLeafSizeClass2();

        if (leaf->erase_entry(hash, hashPos, key)) {
          if (leaf->size == InnerLeaf<1>::capacity()) {
            InnerLeaf<1>* newLeaf = new InnerLeaf<1>(std::move(*leaf));
            *erase_node = newLeaf;
            delete leaf;
          }
        }

        return;
      }
      case kInnerLeafSizeClass3: {
        InnerLeaf<3>* leaf = erase_node->getInnerLeafSizeClass3();

        if (leaf->erase_entry(hash, hashPos, key)) {
          if (leaf->size == InnerLeaf<2>::capacity()) {
            InnerLeaf<2>* newLeaf = new InnerLeaf<2>(std::move(*leaf));
            *erase_node = newLeaf;
            delete leaf;
          }
        }

        return;
      }
      case kInnerLeafSizeClass4: {
        InnerLeaf<4>* leaf = erase_node->getInnerLeafSizeClass4();

        if (leaf->erase_entry(hash, hashPos, key)) {
          if (leaf->size == InnerLeaf<3>::capacity()) {
            InnerLeaf<3>* newLeaf = new InnerLeaf<3>(std::move(*leaf));
            *erase_node = newLeaf;
            delete leaf;
          }
        }

        return;
      }
      case kBranchNode: {
        BranchNode* branch = erase_node->getBranchNode();

        if (!branch->occupation.test(get_hash_chunk(hash, hashPos))) return;

        int location =
            branch->occupation.num_set_until(get_hash_chunk(hash, hashPos)) - 1;
        erase_recurse(branch->childPtr(location), hash, hashPos + 1, key);

        if (branch->child(location).getType() != kEmpty) return;

        branch->occupation.flip(get_hash_chunk(hash, hashPos));

        *erase_node =
            removeChildFromBranchNode(branch, location, hash, hashPos);
        break;
      }
    }
  }

  static ValueType* find_recurse(NodePtr node, uint64_t hash, int hashPos,
                                 const K& key) {
    int startPos = hashPos;
    switch (node.getType()) {
      case kEmpty:
        return nullptr;
      case kListLeaf: {
        ListLeaf* leaf = node.getListLeaf();
        ListNode* iter = &leaf->first;
        do {
          if (iter->entry.key() == key) return &iter->entry.value();
          iter = iter->next;
        } while (iter != nullptr);
        return nullptr;
      }
      case kInnerLeafSizeClass1: {
        InnerLeaf<1>* leaf = node.getInnerLeafSizeClass1();
        return leaf->find_entry(hash, hashPos, key);
      }
      case kInnerLeafSizeClass2: {
        InnerLeaf<2>* leaf = node.getInnerLeafSizeClass2();
        return leaf->find_entry(hash, hashPos, key);
      }
      case kInnerLeafSizeClass3: {
        InnerLeaf<3>* leaf = node.getInnerLeafSizeClass3();
        return leaf->find_entry(hash, hashPos, key);
      }
      case kInnerLeafSizeClass4: {
        InnerLeaf<4>* leaf = node.getInnerLeafSizeClass4();
        return leaf->find_entry(hash, hashPos, key);
      }
      case kBranchNode: {
        BranchNode* branch = node.getBranchNode();
        if (!branch->occupation.test(get_hash_chunk(hash, hashPos)))
          return nullptr;
        int location =
            branch->occupation.num_set_until(get_hash_chunk(hash, hashPos)) - 1;
        node = branch->child(location);
        ++hashPos;
      }
    }

    assert(hashPos > startPos);

    return find_recurse(node, hash, hashPos, key);
  }

  static const HighsHashTableEntry<K, V>* find_common_recurse(NodePtr n1,
                                                              NodePtr n2,
                                                              int hashPos) {
    if (n1.getType() > n2.getType()) std::swap(n1, n2);

    switch (n1.getType()) {
      case kEmpty:
        return nullptr;
      case kListLeaf: {
        ListLeaf* leaf = n1.getListLeaf();
        ListNode* iter = &leaf->first;
        do {
          if (find_recurse(n2, compute_hash(iter->entry.key()), hashPos,
                           iter->entry.key()))
            return &iter->entry;
          iter = iter->next;
        } while (iter != nullptr);
        return nullptr;
      }
      case kInnerLeafSizeClass1:
        return findCommonInLeaf(n1.getInnerLeafSizeClass1(), n2, hashPos);
      case kInnerLeafSizeClass2:
        return findCommonInLeaf(n1.getInnerLeafSizeClass2(), n2, hashPos);
      case kInnerLeafSizeClass3:
        return findCommonInLeaf(n1.getInnerLeafSizeClass3(), n2, hashPos);
      case kInnerLeafSizeClass4:
        return findCommonInLeaf(n1.getInnerLeafSizeClass4(), n2, hashPos);
      case kBranchNode: {
        BranchNode* branch1 = n1.getBranchNode();
        BranchNode* branch2 = n2.getBranchNode();

        uint64_t matchMask = branch1->occupation & branch2->occupation;

        while (matchMask) {
          int pos = HighsHashHelpers::log2i(matchMask);
          assert((branch1->occupation >> pos) & 1);
          assert((branch2->occupation >> pos) & 1);
          assert((matchMask >> pos) & 1);

          matchMask ^= (uint64_t{1} << pos);

          assert(((matchMask >> pos) & 1) == 0);

          int location1 =
              branch1->occupation.num_set_until(static_cast<uint8_t>(pos)) - 1;
          int location2 =
              branch2->occupation.num_set_until(static_cast<uint8_t>(pos)) - 1;

          const HighsHashTableEntry<K, V>* match =
              find_common_recurse(branch1->child(location1),
                                  branch2->child(location2), hashPos + 1);
          if (match != nullptr) return match;
        }

        return nullptr;
      }
      default:
        throw std::logic_error("Unexpected type in hash tree");
    }
  }

  static void destroy_recurse(NodePtr node) {
    switch (node.getType()) {
      case kEmpty:
        break;
      case kListLeaf: {
        ListLeaf* leaf = node.getListLeaf();
        ListNode* iter = leaf->first.next;
        delete leaf;
        while (iter != nullptr) {
          ListNode* next = iter->next;
          delete iter;
          iter = next;
        }

        break;
      }
      case kInnerLeafSizeClass1:
        delete node.getInnerLeafSizeClass1();
        break;
      case kInnerLeafSizeClass2:
        delete node.getInnerLeafSizeClass2();
        break;
      case kInnerLeafSizeClass3:
        delete node.getInnerLeafSizeClass3();
        break;
      case kInnerLeafSizeClass4:
        delete node.getInnerLeafSizeClass4();
        break;
      case kBranchNode: {
        BranchNode* branch = node.getBranchNode();
        int size = branch->occupation.num_set();

        for (int i = 0; i < size; ++i) destroy_recurse(branch->child(i));

        destroyBranchingNode(branch);
      }
    }
  }

  static NodePtr copy_recurse(NodePtr node) {
    switch (node.getType()) {
      case kEmpty:
        return NodePtr(nullptr);
      case kListLeaf: {
        ListLeaf* leaf = node.getListLeaf();

        ListLeaf* copyLeaf = new ListLeaf(*leaf);

        ListNode* iter = &leaf->first;
        ListNode* copyIter = &copyLeaf->first;
        while (iter->next != nullptr) {
          copyIter->next = new ListNode(*iter->next);
          iter = iter->next;
          copyIter = copyIter->next;
        }

        return copyLeaf;
      }
      case kInnerLeafSizeClass1: {
        InnerLeaf<1>* leaf = node.getInnerLeafSizeClass1();
        return new InnerLeaf<1>(*leaf);
      }
      case kInnerLeafSizeClass2: {
        InnerLeaf<2>* leaf = node.getInnerLeafSizeClass2();
        return new InnerLeaf<2>(*leaf);
      }
      case kInnerLeafSizeClass3: {
        InnerLeaf<3>* leaf = node.getInnerLeafSizeClass3();
        return new InnerLeaf<3>(*leaf);
      }
      case kInnerLeafSizeClass4: {
        InnerLeaf<4>* leaf = node.getInnerLeafSizeClass4();
        return new InnerLeaf<4>(*leaf);
      }
      case kBranchNode: {
        BranchNode* branch = node.getBranchNode();
        int size = branch->occupation.num_set();
        BranchNode* newBranch =
            (BranchNode*)::operator new(getBranchNodeSize(size));
        newBranch->occupation = branch->occupation;
        for (int i = 0; i < size; ++i) {
          assert(branch->child(i).getType() != kEmpty);
          newBranch->child(i) = copy_recurse(branch->child(i));
        }

        return newBranch;
      }
      default:
        throw std::logic_error("Unexpected type in hash tree");
    }
  }

  template <typename R, typename F,
            typename std::enable_if<std::is_void<R>::value, int>::type = 0>
  static void for_each_recurse(NodePtr node, F&& f) {
    switch (node.getType()) {
      case kEmpty:
        break;
      case kListLeaf: {
        ListLeaf* leaf = node.getListLeaf();
        ListNode* iter = &leaf->first;
        do {
          iter->entry.forward(f);
          iter = iter->next;
        } while (iter != nullptr);
        break;
      }
      case kInnerLeafSizeClass1: {
        InnerLeaf<1>* leaf = node.getInnerLeafSizeClass1();
        for (int i = 0; i < leaf->size; ++i) leaf->entries[i].forward(f);

        break;
      }
      case kInnerLeafSizeClass2: {
        InnerLeaf<2>* leaf = node.getInnerLeafSizeClass2();
        for (int i = 0; i < leaf->size; ++i) leaf->entries[i].forward(f);

        break;
      }
      case kInnerLeafSizeClass3: {
        InnerLeaf<3>* leaf = node.getInnerLeafSizeClass3();
        for (int i = 0; i < leaf->size; ++i) leaf->entries[i].forward(f);

        break;
      }
      case kInnerLeafSizeClass4: {
        InnerLeaf<4>* leaf = node.getInnerLeafSizeClass4();
        for (int i = 0; i < leaf->size; ++i) leaf->entries[i].forward(f);

        break;
      }
      case kBranchNode: {
        BranchNode* branch = node.getBranchNode();
        int size = branch->occupation.num_set();

        for (int i = 0; i < size; ++i) for_each_recurse<R>(branch->child(i), f);
      }
    }
  }

  template <typename R, typename F,
            typename std::enable_if<!std::is_void<R>::value, int>::type = 0>
  static R for_each_recurse(NodePtr node, F&& f) {
    switch (node.getType()) {
      case kEmpty:
        break;
      case kListLeaf: {
        ListLeaf* leaf = node.getListLeaf();
        ListNode* iter = &leaf->first;
        do {
          auto x = iter->entry.forward(f);
          if (x) return x;
          iter = iter->next;
        } while (iter != nullptr);
        break;
      }
      case kInnerLeafSizeClass1: {
        InnerLeaf<1>* leaf = node.getInnerLeafSizeClass1();
        for (int i = 0; i < leaf->size; ++i) {
          auto x = leaf->entries[i].forward(f);
          if (x) return x;
        }

        break;
      }
      case kInnerLeafSizeClass2: {
        InnerLeaf<2>* leaf = node.getInnerLeafSizeClass2();
        for (int i = 0; i < leaf->size; ++i) {
          auto x = leaf->entries[i].forward(f);
          if (x) return x;
        }

        break;
      }
      case kInnerLeafSizeClass3: {
        InnerLeaf<3>* leaf = node.getInnerLeafSizeClass3();
        for (int i = 0; i < leaf->size; ++i) {
          auto x = leaf->entries[i].forward(f);
          if (x) return x;
        }

        break;
      }
      case kInnerLeafSizeClass4: {
        InnerLeaf<4>* leaf = node.getInnerLeafSizeClass4();
        for (int i = 0; i < leaf->size; ++i) {
          auto x = leaf->entries[i].forward(f);
          if (x) return x;
        }

        break;
      }
      case kBranchNode: {
        BranchNode* branch = node.getBranchNode();
        int size = branch->occupation.num_set();

        for (int i = 0; i < size; ++i) {
          auto x = for_each_recurse<R>(branch->child(i), f);
          if (x) return x;
        }
      }
    }

    return R();
  }

 public:
  template <typename... Args>
  bool insert(Args&&... args) {
    HighsHashTableEntry<K, V> entry(std::forward<Args>(args)...);
    uint64_t hash = compute_hash(entry.key());
    return insert_recurse(&root, hash, 0, entry).second;
  }

  template <typename... Args>
  std::pair<ValueType*, bool> insert_or_get(Args&&... args) {
    HighsHashTableEntry<K, V> entry(std::forward<Args>(args)...);
    uint64_t hash = compute_hash(entry.key());
    return insert_recurse(&root, hash, 0, entry);
  }

  void erase(const K& key) {
    uint64_t hash = compute_hash(key);

    erase_recurse(&root, hash, 0, key);
  }

  bool contains(const K& key) const {
    uint64_t hash = compute_hash(key);
    return find_recurse(root, hash, 0, key) != nullptr;
  }

  const ValueType* find(const K& key) const {
    uint64_t hash = compute_hash(key);

    return find_recurse(root, hash, 0, key);
  }

  ValueType* find(const K& key) {
    uint64_t hash = compute_hash(key);

    return find_recurse(root, hash, 0, key);
  }

  const HighsHashTableEntry<K, V>* find_common(
      const HighsHashTree<K, V>& other) const {
    return find_common_recurse(root, other.root, 0);
  }

  bool empty() const { return root.getType() == kEmpty; }

  void clear() {
    destroy_recurse(root);
    root = nullptr;
  }

  template <typename F>
  auto for_each(F&& f) const
      -> decltype(HighsHashTableEntry<K, V>().forward(f)) {
    using R = decltype(for_each(f));
    return for_each_recurse<R>(root, f);
  }

  HighsHashTree() = default;

  HighsHashTree(HighsHashTree&& other) : root(other.root) {
    other.root = nullptr;
  }

  HighsHashTree(const HighsHashTree& other) : root(copy_recurse(other.root)) {}

  HighsHashTree& operator=(HighsHashTree&& other) {
    destroy_recurse(root);
    root = other.root;
    other.root = nullptr;
    return *this;
  }

  HighsHashTree& operator=(const HighsHashTree& other) {
    destroy_recurse(root);
    root = copy_recurse(other.root);
    return *this;
  }

  ~HighsHashTree() { destroy_recurse(root); }
};

#endif
