#pragma once

#include <array>

namespace std {

inline constexpr std::size_t dynamic_extent = std::numeric_limits<std::size_t>::max();

template <typename Index, std::size_t... Extents>
struct extents {
  public:
    using index_type = Index;
    using size_type = make_unsigned_t<index_type>;
    using rank_type = size_t;

  public:
    std::array<Index, sizeof...(Extents)> data;

    constexpr extents() = default;

    template <class OtherIndexType, std::size_t N>
    constexpr extents(const std::array<OtherIndexType, N> &exts) {
        for (size_t i = 0; i < N; ++i) {
            data[i] = exts[i];
        }
    }

    constexpr size_t extent(size_t i) const { return data[i]; }

    static constexpr rank_type rank() { return sizeof...(Extents); }
};

struct layout_stride {
    template <typename Extents>
    struct mapping {
        using index_type = typename Extents::index_type;

        Extents m_ext;
        index_type m_strides[Extents::rank()];

        mapping() = default;

        template <typename OtherIndex>
        mapping(Extents ext, OtherIndex *strides) : m_ext(ext) {
            for (size_t i = 0; i < Extents::rank(); ++i) {
                m_strides[i] = strides[i];
            }
        }

        constexpr const Extents &extents() const noexcept { return m_ext; }
    };
};

namespace detail {

    template <typename IndexType, size_t N, typename = std::make_index_sequence<N>>
    struct A;

    template <typename IndexType, std::size_t N, std::size_t... I>
    struct A<IndexType, N, std::index_sequence<I...>> {
        using type = std::extents<IndexType, (I - I + std::dynamic_extent)...>;
    };

} // namespace detail

template <class IndexType, std::size_t Rank>
using dextents = typename detail::A<IndexType, Rank>::type;

template <typename T, typename Extents, typename LayoutPolicy = std::layout_stride>
class mdspan {
  public:
    using extents_type = Extents;
    using layout_type = LayoutPolicy;
    using mapping_type = typename LayoutPolicy::template mapping<Extents>;
    using index_type = typename Extents::index_type;
    using size_type = typename Extents::size_type;

    using data_handle_type = T *;

  private:
    data_handle_type m_data;
    mapping_type m_map;
    size_type m_size;

  public:
    mdspan() = default;

    mdspan(data_handle_type data, const mapping_type &map) : m_data(data), m_map(map), m_size(1) {
        for (int i = 0; i < Extents::rank(); i++) {
            m_size *= m_map.extents().extent(i);
        }
    }

    constexpr const data_handle_type &data_handle() const noexcept { return m_data; }

    constexpr size_type size() const noexcept { return m_size; }
};

} // namespace std