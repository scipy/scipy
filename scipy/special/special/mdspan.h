#pragma once

#include <array>
#include <limits>
#include <tuple>
#include <type_traits>
#include <utility>

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

    template <class... OtherIndex>
    constexpr explicit extents(OtherIndex... exts) : m_dexts{exts...} {}

    template <class OtherIndex>
    constexpr extents(const array<OtherIndex, sizeof...(Extents)> &dexts) noexcept : m_dexts(dexts) {}

    constexpr index_type extent(rank_type i) const noexcept { return m_dexts[i]; }

    static constexpr rank_type rank() noexcept { return sizeof...(Extents); }
};

template <typename Index, size_t Rank>
using dextents = detail::fill_extents_t<Index, Rank, dynamic_extent>;

struct full_extent_t {
    explicit full_extent_t() = default;
};

inline constexpr full_extent_t full_extent;

template <typename Offset, typename Extent, typename Stride>
struct strided_slice {
    using offset_type = Offset;
    using extent_type = Extent;
    using stride_type = Stride;

    strided_slice() = default;

    strided_slice(offset_type offset, extent_type extent, stride_type stride)
        : offset(offset), extent(extent), stride(stride) {}

    offset_type offset;
    extent_type extent;
    stride_type stride;
};

namespace detail {

    template <typename Index, typename Offset, typename Extent, typename Stride>
    Index submdspan_extent(Index ext, strided_slice<Offset, Extent, Stride> slice) {
        return (slice.extent - slice.offset) / slice.stride;
    }

    template <typename Index, typename Index0, typename Index1>
    Index submdspan_extent(Index ext, std::tuple<Index0, Index1> slice) {
        return std::get<1>(slice) - std::get<0>(slice);
    }

    template <typename Index>
    Index submdspan_extent(Index ext, full_extent_t slice) {
        return ext;
    }

    template <size_t... I, typename Index, size_t... Extents, typename... Slices>
    auto submdspan_extents(std::index_sequence<I...>, const extents<Index, Extents...> exts, Slices... slices) {
        return extents<Index, Extents...>{submdspan_extent(exts.extent(I), slices)...};
    }

    template <typename Index, size_t... Extents, typename... Slices>
    auto submdspan_extents(const extents<Index, Extents...> exts, Slices... slices) {
        return submdspan_extents(std::index_sequence_for<Slices...>(), exts, slices...);
    }

} // namespace detail

template <typename Index, size_t... Extents, typename... Slices>
auto submdspan_extents(const extents<Index, Extents...> &exts, Slices... slices) {
    return detail::submdspan_extents(exts, slices...);
}

template <typename Mapping, typename... Slices>
auto submdspan_mapping(const Mapping &, Slices...);

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

        constexpr const extents_type &extents() const noexcept { return m_exts; }

        constexpr const array<index_type, extents_type::rank()> &strides() const noexcept { return m_strides; }

        constexpr index_type extent(rank_type i) const noexcept { return m_exts.extent(i); }

        constexpr index_type stride(rank_type i) const noexcept { return m_strides[i]; }

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

namespace detail {

    template <typename Index, typename Offset, typename Extent, typename Stride>
    Index submdspan_stride(Index stride, strided_slice<Offset, Extent, Stride> slice) {
        return stride * slice.stride;
    }

    template <typename Index, typename Index0, typename Index1>
    Index submdspan_stride(Index stride, std::tuple<Index0, Index1> slice) {
        return stride;
    }

    template <typename Index>
    Index submdspan_stride(Index stride, full_extent_t slice) {
        return stride;
    }

    template <size_t... I, typename Index, size_t Rank, typename... Slices>
    auto submdspan_strides(std::index_sequence<I...>, const array<Index, Rank> strides, Slices... slices) {
        array<Index, Rank> res{submdspan_stride(strides[I], slices)...};
        return res;
    }

    template <typename Index, size_t Rank, typename... Slices>
    auto submdspan_strides(const array<Index, Rank> strides, Slices... slices) {
        return submdspan_strides(std::index_sequence_for<Slices...>(), strides, slices...);
    }

} // namespace detail

template <typename Index, size_t Rank, typename... Slices>
auto submdspan_strides(const array<Index, Rank> &strides, Slices... slices) {
    return detail::submdspan_strides(strides, slices...);
}

template <typename Extents, typename... Slices>
auto submdspan_mapping(const layout_stride::mapping<Extents> &map, Slices... slices) {
    return layout_stride::mapping(submdspan_extents(map.extents(), slices...),
                                  submdspan_strides(map.strides(), slices...));
}

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

    constexpr mdspan(data_handle_type p, const mapping_type &m, const accessor_type &a)
        : m_ptr(p), m_map(m), m_acc(a) {}

    template <typename... OtherIndices>
    constexpr reference operator()(OtherIndices... indices) const {
        return m_acc.access(m_ptr, m_map(static_cast<index_type>(std::move(indices))...));
    }

    template <typename OtherIndex>
    constexpr reference operator[](OtherIndex index) const {
        return m_acc.access(m_ptr, m_map(static_cast<index_type>(index)));
    }

    constexpr const data_handle_type &data_handle() const noexcept { return m_ptr; }

    constexpr const mapping_type &mapping() const noexcept { return m_map; }

    constexpr const accessor_type &accessor() const noexcept { return m_acc; }

    constexpr index_type stride(rank_type r) const { return m_map.stride(r); }

    constexpr const extents_type &extents() const noexcept { return m_map.extents(); }

    constexpr index_type extent(rank_type r) const noexcept { return m_map.extent(r); }

    constexpr size_type size() const noexcept {
        size_type res = 1;
        for (rank_type i = 0; i < extents_type::rank(); ++i) {
            res *= m_map.extent(i);
        }

        return res;
    }
};

namespace detail {

    template <typename Offset, typename Extent, typename Stride>
    auto submdspan_offset(strided_slice<Offset, Extent, Stride> slice) {
        return slice.offset;
    }

    template <typename Index0, typename Index1>
    auto submdspan_offset(std::tuple<Index0, Index1> slice) {
        return std::get<0>(slice);
    }

    inline auto submdspan_offset(full_extent_t slice) { return 0; }

} // namespace detail

template <class T, class Extents, class LayoutPolicy, class AccessorPolicy, class... SliceArgs>
auto submdspan(const mdspan<T, Extents, LayoutPolicy, AccessorPolicy> &src, SliceArgs... args) {
    static_assert(Extents::rank() == sizeof...(SliceArgs), "number of slices must equal extents rank");

    using submdspan_type = mdspan<T, Extents, LayoutPolicy, AccessorPolicy>;

    auto src_map = src.mapping();
    auto src_acc = src.accessor();
    return submdspan_type(src_acc.offset(src.data_handle(), src_map(detail::submdspan_offset(args)...)),
                          submdspan_mapping(src.mapping(), args...), src_acc);
}

} // namespace std
