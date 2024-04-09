#pragma once

#include <array>
#include <limits>
#include <type_traits>

namespace std {

template <typename Index, size_t... Extents>
class extents;

namespace detail {

    template <typename Index, size_t Rank, size_t Value, size_t... Extents>
    struct fill_extents {
        using type = typename fill_extents<Index, Rank - 1, Value, Extents..., Value>::type;
    };

    template <typename Index, size_t Value, size_t... Extents>
    struct fill_extents<Index, 0, Value, Extents...> {
        using type = extents<Index, Extents...>;
    };

    template <typename Index, size_t Rank, size_t Value>
    using fill_extents_t = typename fill_extents<Index, Rank, Value>::type;

} // namespace detail

inline constexpr size_t dynamic_extent = numeric_limits<size_t>::max();

template <typename Index, size_t... Extents>
class extents {
    static_assert(((Extents == dynamic_extent) && ... && true), "extents must all be dynamic");

  public:
    using index_type = Index;
    using size_type = make_unsigned_t<index_type>;
    using rank_type = size_t;

  private:
    array<index_type, sizeof...(Extents)> m_dexts;

  public:
    constexpr extents() = default;

    template <class OtherIndex>
    constexpr extents(const array<OtherIndex, sizeof...(Extents)> &dexts) noexcept : m_dexts(dexts) {}

    constexpr index_type extent(rank_type i) const noexcept { return m_dexts[i]; }

    static constexpr rank_type rank() noexcept { return sizeof...(Extents); }
};

template <typename Index, size_t Rank>
using dextents = detail::fill_extents_t<Index, Rank, dynamic_extent>;

struct layout_left;

struct layout_right;

struct layout_stride {
    template <typename Extents>
    class mapping {
      public:
        using extents_type = Extents;
        using index_type = typename extents_type::index_type;
        using size_type = typename extents_type::size_type;
        using rank_type = typename extents_type::rank_type;
        using layout_type = layout_stride;

      private:
        extents_type m_exts;
        array<index_type, extents_type::rank()> m_strides;

      public:
        constexpr mapping() = default;

        constexpr mapping(const Extents &exts, const array<index_type, extents_type::rank()> &strides)
            : m_exts(exts), m_strides(strides) {}

        constexpr index_type extent(rank_type i) const noexcept { return m_exts.extent(i); }

        template <typename... Args>
        constexpr index_type operator()(Args... args) const noexcept {
            static_assert(sizeof...(Args) == extents_type::rank(), "index must have same rank as extents");

            index_type indices[extents_type::rank()] = {args...};
            index_type res = 0;
            for (rank_type i = 0; i < extents_type::rank(); ++i) {
                res += indices[i] * m_strides[i];
            }

            return res;
        }
    };
};

template <typename Element>
class default_accessor {
  public:
    using offset_policy = default_accessor;
    using element_type = Element;
    using reference = Element &;
    using data_handle_type = Element *;

    constexpr reference access(data_handle_type p, size_t i) const noexcept { return p[i]; }

    constexpr data_handle_type offset(data_handle_type p, size_t i) const noexcept { return p + i; }
};

template <typename T, typename Extents, typename LayoutPolicy = layout_right,
          typename AccessorPolicy = default_accessor<T>>
class mdspan {
  public:
    using extents_type = Extents;
    using layout_type = LayoutPolicy;
    using accessor_type = AccessorPolicy;
    using mapping_type = typename LayoutPolicy::template mapping<Extents>;
    using element_type = T;
    using value_type = remove_cv_t<T>;
    using index_type = typename Extents::index_type;
    using size_type = typename Extents::size_type;
    using rank_type = typename Extents::rank_type;
    using data_handle_type = typename AccessorPolicy::data_handle_type;
    using reference = typename AccessorPolicy::reference;

  private:
    data_handle_type m_ptr;
    mapping_type m_map;
    accessor_type m_acc;

  public:
    constexpr mdspan() = default;

    constexpr mdspan(data_handle_type p, const mapping_type &m) : m_ptr(p), m_map(m) {}

    template <typename... OtherIndices>
    constexpr reference operator()(OtherIndices... indices) const {
        return m_acc.access(m_ptr, m_map(static_cast<index_type>(std::move(indices))...));
    }

    template <typename OtherIndex>
    constexpr reference operator[](OtherIndex index) const {
        return m_acc.access(m_ptr, m_map(static_cast<index_type>(index)));
    }

    constexpr index_type extent(rank_type i) const noexcept { return m_map.extent(i); }

    constexpr size_type size() const noexcept {
        size_type res = 1;
        for (rank_type i = 0; i < extents_type::rank(); ++i) {
            res *= m_map.extent(i);
        }

        return res;
    }
};

} // namespace std
