
#pragma once

namespace std {

inline constexpr std::size_t dynamic_extent = std::numeric_limits<std::size_t>::max();

template <typename Index, size_t... Extents>
class extents {
  public:
    using index_type = Index;
    using size_type = make_unsigned_t<index_type>;
    using rank_type = size_t;

  public:
    extents() = default;

    static constexpr rank_type rank() noexcept { return sizeof...(Extents); }

    static constexpr rank_type rank_dynamic() noexcept { return rank(); }

    template <class OtherIndexType>
    constexpr extents(const std::array<OtherIndexType, rank_dynamic()> &dexts) noexcept : m_dexts(dexts) {}

    constexpr index_type extent(rank_type i) const noexcept { return m_dexts[i]; }

    constexpr size_type size() const noexcept {
        size_type res = 1;
        for (size_t i = 0; i < m_dexts.size(); ++i) {
            res *= m_dexts[i];
        }

        return res;
    }

  private:
    std::array<index_type, rank_dynamic()> m_dexts;
};

template <typename IndexType, std::size_t Rank>
using dextents = std::conditional_t<Rank == 1, extents<IndexType, dynamic_extent>,
                                    extents<IndexType, dynamic_extent, dynamic_extent>>;

struct layout_right;

struct layout_stride {
    template <typename Extents>
    struct mapping {
      public:
        using extents_type = Extents;
        using index_type = typename extents_type::index_type;
        using size_type = typename extents_type::size_type;
        using rank_type = typename extents_type::rank_type;

      public:
        extents_type m_exts;
        array<index_type, extents_type::rank()> m_strs;

      public:
        mapping() = default;

        template <typename OtherExtents>
        mapping(const OtherExtents &exts, const array<index_type, extents_type::rank()> &s) : m_exts(exts), m_strs(s) {}

        template <typename... Index>
        constexpr index_type operator()(Index... indices) const noexcept {
            std::array<index_type, sizeof...(Index)> ind{indices...};
            return operator()(ind);
        }

        template <class OtherIndex>
        constexpr index_type operator()(const std::array<OtherIndex, extents_type::rank()> &indices) const {
            index_type res = 0;
            for (size_t i = 0; i < indices.size(); ++i) {
                res += indices[i] * m_strs[i];
            }

            return res;
        }

        constexpr index_type extent(rank_type i) const noexcept { return m_exts.extent(i); }

        constexpr size_type size() const noexcept { return m_exts.size(); }
    };
};

template <typename ElementType>
struct default_accessor {
    using reference = ElementType &;
    using data_handle_type = ElementType *;

    constexpr reference access(data_handle_type p, std::size_t i) const noexcept { return p[i]; }
};

template <class T, class Extents, class LayoutPolicy = std::layout_right,
          class AccessorPolicy = std::default_accessor<T>>
class mdspan {
  public:
    using extents_type = Extents;
    using accessor_type = AccessorPolicy;
    using mapping_type = typename LayoutPolicy::template mapping<Extents>;
    using reference = typename AccessorPolicy::reference;
    using index_type = typename Extents::index_type;
    using size_type = typename Extents::size_type;
    using rank_type = typename Extents::rank_type;
    using data_handle_type = typename AccessorPolicy::data_handle_type;

  private:
    data_handle_type m_ptr;
    mapping_type m_map;
    accessor_type m_acc;

  public:
    constexpr mdspan() = default;

    constexpr mdspan(data_handle_type p, const mapping_type &m) : m_ptr(p), m_map(m) {}

    template <class... OtherIndexTypes>
    constexpr reference operator()(OtherIndexTypes... indices) const {
        return m_acc.access(m_ptr, m_map(static_cast<index_type>(std::move(indices))...));
    }

    template <class OtherIndexType>
    constexpr reference operator[](OtherIndexType index) const {
        return m_acc.access(m_ptr, m_map(static_cast<index_type>(index)));
    }

    constexpr index_type extent(rank_type i) const noexcept { return m_map.extent(i); }

    constexpr size_type size() const noexcept { return m_map.size(); }
};

} // namespace std
