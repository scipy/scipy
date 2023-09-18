/****************************************************************************
**
** Copyright (C) 2014 Klaralvdalens Datakonsult AB (KDAB).
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt3D module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:LGPL$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** GNU Lesser General Public License Usage
** Alternatively, this file may be used under the terms of the GNU Lesser
** General Public License version 3 as published by the Free Software
** Foundation and appearing in the file LICENSE.LGPL3 included in the
** packaging of this file. Please review the following information to
** ensure the GNU Lesser General Public License version 3 requirements
** will be met: https://www.gnu.org/licenses/lgpl-3.0.html.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 2.0 or (at your option) the GNU General
** Public license version 3 or any later version approved by the KDE Free
** Qt Foundation. The licenses are as published by the Free Software
** Foundation and appearing in the file LICENSE.GPL2 and LICENSE.GPL3
** included in the packaging of this file. Please review the following
** information to ensure the GNU General Public License requirements will
** be met: https://www.gnu.org/licenses/gpl-2.0.html and
** https://www.gnu.org/licenses/gpl-3.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef QT3DCORE_QCIRCULARBUFFER_H
#define QT3DCORE_QCIRCULARBUFFER_H

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API.  It exists for the convenience
// of other Qt classes.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.
//

#include <Qt3DCore/qt3dcore_global.h>
#include <QtCore/QtGlobal>
#include <QtCore/qlist.h>
#include <QtCore/qpair.h>
#include <QtCore/qshareddata.h>
#include <QtCore/qvector.h>

#include <algorithm>
#include <iterator>
#include <limits>
#include <memory>
#include <new>


#ifdef Q_COMPILER_INITIALIZER_LISTS
# include <initializer_list>
#endif

QT_BEGIN_NAMESPACE

namespace Qt3DCore {

class CircularBufferData : public QSharedData
{
protected:
    CircularBufferData()
        : data(nullptr),
          capacity(0),
          size(0),
          first(-1),
          last(-1)
    {}

    ~CircularBufferData()
    {
        // Release the raw memory
        deallocate(data);
    }

    int wraparound(int index) const
    {
        index = index < capacity ? index : (index - capacity);
        Q_ASSERT(index < capacity); // make sure that wrapping around once was enough
        return index;
    }

    void *allocate(int count, size_t sizeOfT)
    { return operator new[](count * sizeOfT); }
    void deallocate(void *p)
    { operator delete[](p); }

    void *data;      // Array of the actual data
public:
    int capacity; // Size of the m_data array
    int size;     // Number of elements of m_data actually used
    int first;    // Index in m_data of the first element of the circular buffer
    int last;     // Index in m_data of the last element of the circular buffer
};

template <typename T>
class TypedCircularBufferData : public CircularBufferData
{
    template <typename InputIterator>
    explicit TypedCircularBufferData(InputIterator f, InputIterator l, std::input_iterator_tag) Q_DECL_EQ_DELETE;
public:
    TypedCircularBufferData() : CircularBufferData() {}
    template <typename ForwardIterator>
    explicit TypedCircularBufferData(ForwardIterator f, ForwardIterator l, std::forward_iterator_tag)
    {
        const int n = int(std::distance(f, l));
        CircularBufferData::data = allocate(n);
        std::uninitialized_copy(f, l, data());
        first = 0;
        last = n - 1;
        size = capacity = n;
    }
    ~TypedCircularBufferData()
    {
        if (QTypeInfo<T>::isComplex && size != 0) {
            // The type is complex so we manually call the destructor for each item
            // since we used the placement new operator to instantiate them
            if (first <= last) {
                // Destroy items from first to last
                T *b = data() + first;
                T *i = b + size;
                while (i-- != b)
                     i->~T();
            } else {
                // Destroy items at end of data array
                T *b = data() + first;
                T *i = data() + capacity;
                while (i-- != b)
                     i->~T();

                // Destroy items at beginning of data array
                b = data();
                i = b + last;
                while (i-- != b)
                     i->~T();
            }
        }

    }

    using CircularBufferData::wraparound;
    T *allocate(int count) { return static_cast<T*>(CircularBufferData::allocate(count, sizeof(T))); }
    using CircularBufferData::deallocate;
    T *data() const { return static_cast<T*>(CircularBufferData::data); }
    void setData(T *newdata) { CircularBufferData::data = static_cast<void*>(newdata); }
};

template <typename T>
class QCircularBuffer
{
    typedef TypedCircularBufferData<T> Data;
public:
    typedef QPair<T*,int> array_range;
    typedef QPair<const T*,int> const_array_range;
    typedef array_range ArrayRange;
    typedef const_array_range ConstArrayRange;

    QCircularBuffer()
        : d(new Data())
    {}

    explicit QCircularBuffer(int amount);
    explicit QCircularBuffer(int amount, const T &val);
    explicit QCircularBuffer(int amount, int initialSize, const T &value);
#ifdef Q_COMPILER_INITIALIZER_LISTS
    QCircularBuffer(std::initializer_list<T> list)
        : d(new Data(list.begin(), list.end(), std::random_access_iterator_tag()))
    {}
#endif
    template <typename ForwardIterator>
    explicit QCircularBuffer(ForwardIterator f, ForwardIterator l)
        : d(new Data(f, l, typename std::iterator_traits<ForwardIterator>::iterator_category()))
    {}

    QCircularBuffer(const QCircularBuffer<T> &other)
        : d(other.d)
    {}

    void swap(QCircularBuffer &other) { d.swap(other.d); }

    QCircularBuffer<T> &operator = (const QCircularBuffer<T> &other)
    {
        d = other.d;
        return *this;
    }

    ~QCircularBuffer() {}

    class iterator
    {
    public:
        typedef std::random_access_iterator_tag iterator_category;
        typedef ptrdiff_t difference_type;
        typedef T value_type;
        typedef T *pointer;
        typedef T &reference;

        Q_DECL_CONSTEXPR iterator() : buffer(nullptr), index(-1) {}
        iterator(QCircularBuffer<T> *buf, int idx)
            : buffer(buf), index(idx)
        {}

        T &operator*() const { return (*buffer)[ index ]; }
        T *operator->() const { return &(*buffer)[index]; }
        T &operator[](int j) const { return (*buffer)[ index + j ]; }

        bool operator==(const iterator &other) const
        {
            return (buffer == other.buffer && index == other.index);
        }

        bool operator!=(const iterator &other) const
        {
            return (buffer != other.buffer || index != other.index);
        }

        bool operator<(const iterator &other) const
        {
            Q_ASSERT_X(buffer == other.buffer, "QCircularBuffer<T>::iterator::operator<", "iterators use different buffers");
            return index < other.index;
        }

        bool operator<=(const iterator &other) const
        {
            Q_ASSERT_X(buffer == other.buffer, "QCircularBuffer<T>::iterator::operator<=", "iterators use different buffers");
            return index <= other.index;
        }

        bool operator>(const iterator &other) const
        {
            Q_ASSERT_X(buffer == other.buffer, "QCircularBuffer<T>::iterator::operator>", "iterators use different buffers");
            return index > other.index;
        }

        bool operator>=(const iterator &other) const
        {
            Q_ASSERT_X(buffer == other.buffer, "QCircularBuffer<T>::iterator::operator>=", "iterators use different buffers");
            return index >= other.index;
        }

        iterator &operator++() { ++index; return *this; }
        iterator operator++(int)
        {
            iterator ans = *this;
            ++index;
            return ans;
        }

        iterator &operator--() { --index; return *this; }
        iterator operator--(int)
        {
            iterator ans = *this;
            --index;
            return ans;
        }

        iterator &operator+=(int j) { index += j; return *this; }
        iterator &operator-=(int j) { index -= j; return *this; }
        iterator operator+(int j) const { return iterator(buffer, index + j); }
        iterator operator-(int j) const { return iterator(buffer, index - j); }
        int operator-(iterator other) const
        {
            Q_ASSERT_X(buffer == other.buffer, "QCircularBuffer<T>::iterator::operator-", "iterators use different buffers");
            return index - other.index;
        }

    private:
        QCircularBuffer<T> *buffer;
        int index;
        friend class QCircularBuffer;
    };
    friend class iterator;

    class const_iterator
    {
    public:
        typedef std::random_access_iterator_tag  iterator_category;
        typedef ptrdiff_t difference_type;
        typedef T value_type;
        typedef const T *pointer;
        typedef const T &reference;

        Q_DECL_CONSTEXPR const_iterator() : buffer(nullptr), index(-1) {}
        const_iterator(const QCircularBuffer<T> *buff, int idx)
            : buffer(buff), index(idx)
        {}
        const_iterator(const iterator &other)
            : buffer(other.buffer), index(other.index)
        {}

        const T &operator*() const { return buffer->at(index); }
        const T *operator->() const { return &buffer->at(index); }
        const T &operator[](int j) const { return buffer->at(index + j); }

        bool operator==(const const_iterator &other) const
        {
            Q_ASSERT_X(buffer == other.buffer, "QCircularBuffer<T>::const_iterator::operator==", "iterators use different buffers");
            return index == other.index;
        }

        bool operator!=(const const_iterator &other) const
        {
            Q_ASSERT_X(buffer == other.buffer, "QCircularBuffer<T>::const_iterator::operator!=", "iterators use different buffers");
            return index != other.index;
        }

        bool operator<(const const_iterator &other) const
        {
            Q_ASSERT_X(buffer == other.buffer, "QCircularBuffer<T>::const_iterator::operator<", "iterators use different buffers");
            return index < other.index;
        }

        bool operator<=(const const_iterator &other) const
        {
            Q_ASSERT_X(buffer == other.buffer, "QCircularBuffer<T>::const_iterator::operator<=", "iterators use different buffers");
            return index <= other.index;
        }

        bool operator>(const const_iterator &other) const
        {
            Q_ASSERT_X(buffer == other.buffer, "QCircularBuffer<T>::const_iterator::operator>", "iterators use different buffers");
            return index > other.index;
        }

        bool operator>=(const const_iterator &other) const
        {
            Q_ASSERT_X(buffer == other.buffer, "QCircularBuffer<T>::const_iterator::operator>=", "iterators use different buffers");
            return index >= other.index;
        }

        const_iterator &operator++() { ++index; return *this; }
        const_iterator operator++(int)
        {
            const_iterator ans = *this;
            ++index;
            return ans;
        }

        const_iterator &operator--() { --index; return *this; }
        const_iterator operator--(int)
        {
            const_iterator ans = *this;
            --index;
            return ans;
        }

        const_iterator &operator+=(int j) { index += j; return *this; }
        const_iterator &operator-=(int j) { index -= j; return *this; }
        const_iterator operator+(int j) const { return const_iterator(buffer, index + j); }
        const_iterator operator-(int j) const { return const_iterator(buffer, index - j); }
        int operator-(const_iterator other) const
        {
            Q_ASSERT_X(buffer == other.buffer, "QCircularBuffer<T>::const_iterator::operator-", "iterators use different buffers");
            return index - other.index;
        }

    private:
        const QCircularBuffer<T> *buffer;
        int index;
        friend class QCircularBuffer;
    };
    friend class const_iterator;

    typedef std::reverse_iterator<iterator> reverse_iterator;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

    iterator begin() { return iterator(this, 0); }
    const_iterator begin() const { return const_iterator(this, 0); }
    const_iterator cbegin() const { return const_iterator(this, 0); }
    reverse_iterator rbegin() { return reverse_iterator(end()); }
    const_reverse_iterator rbegin() const { return const_reverse_iterator(end()); }
    const_reverse_iterator crbegin() const { return const_reverse_iterator(end()); }
    const_iterator constBegin() const { return const_iterator(this, 0); }
    iterator end() { return iterator(this, d->size); }
    const_iterator end() const { return const_iterator(this, d->size); }
    const_iterator cend() const { return const_iterator(this, d->size); }
    reverse_iterator rend() { return reverse_iterator(begin()); }
    const_reverse_iterator rend() const { return const_reverse_iterator(begin()); }
    const_reverse_iterator crend() const { return const_reverse_iterator(begin()); }
    const_iterator constEnd() const { return const_iterator(this, d->size); }
    iterator insert(const_iterator before, int number, const T &val)
    {
        insert(before.index, number, val);
        return iterator(this, before.index);
    }
    iterator insert(const_iterator before, const T &val) { return insert(before, 1, val); }
    iterator erase(const_iterator b, const_iterator e)
    {
        int number = e - b;
        remove(b.index, number);
        return iterator(this, e.index - number);
    }
    iterator erase(const_iterator pos) { return erase(pos, pos + 1); }

    // STL compatibility
    typedef T value_type;
    typedef value_type *pointer;
    typedef const value_type *const_pointer;
    typedef value_type &reference;
    typedef const value_type &const_reference;
    typedef ptrdiff_t difference_type;
    typedef iterator Iterator;
    typedef const_iterator ConstIterator;
    typedef int size_type;

    void push_back(const T &t) { append(t); }
    void push_front(const T &t) { prepend(t); }
    void pop_back() { Q_ASSERT(!isEmpty()); erase(end() - 1); }
    void pop_front() { Q_ASSERT(!isEmpty()); erase(begin()); }
    bool empty() const { return isEmpty(); }
    T &front() { return first(); }
    const T &front() const { return first(); }
    T &back() { return last(); }
    const T &back() const { return last(); }

    QAtomicInt refCount() const { return d->ref; }

    void append(const T &val);

    const T &at(int i) const
    {
        Q_ASSERT_X(i >= 0 && i < d->size, "QCircularBuffer<T>::at", "index out of range");
        return d->data()[d->wraparound(d->first + i)];
    }

    const T &operator[](int i) const
    {
        Q_ASSERT_X(i >= 0 && i < d->size, "QCircularBuffer<T>::operator[]", "index out of range");
        return d->data()[d->wraparound(d->first + i)];
    }

    T &operator[](int i)
    {
        d.detach();
        Q_ASSERT_X(i >= 0 && i < d->size, "QCircularBuffer<T>::operator[]", "index out of range");
        return d->data()[d->wraparound(d->first + i)];
    }

    int capacity() const { return d->capacity; }

    void clear() { *this = QCircularBuffer<T>(d->capacity); }

    bool contains(const T &val) const;
    int count(const T &val) const;
    int count() const { return size(); }

    array_range data()
    {
        d.detach();
        if (d->size == 0)
            return array_range(nullptr, 0);
        if (!isLinearised())
            linearise();
        return array_range(d->data() + d->first, d->last - d->first + 1);
    }
    const_array_range data() const { return constData(); }
    const_array_range constData() const
    {
        if (!isLinearised() || d->size == 0)
            return const_array_range(nullptr, 0);
        return const_array_range(d->data() + d->first, d->last - d->first + 1);
    }

    array_range dataOne()
    {
        d.detach();
        if (d->size == 0)
            return array_range(nullptr, 0);
        if (isLinearised())
            return array_range(d->data() + d->first, d->last - d->first + 1);
        else
            return array_range(d->data() + d->first, d->capacity - d->first);
    }
    const_array_range dataOne() const { return constDataOne(); }
    const_array_range constDataOne() const
    {
        if (d->size == 0)
            return const_array_range(nullptr, 0);
        if (isLinearised())
            return const_array_range(d->data() + d->first, d->last - d->first + 1);
        else
            return const_array_range(d->data() + d->first, d->capacity - d->first);
    }

    array_range dataTwo()
    {
        d.detach();
        if (d->size == 0 || isLinearised())
            return array_range(nullptr, 0);
        return array_range(d->data(), d->last + 1);
    }
    const_array_range dataTwo() const { return constDataTwo(); }
    const_array_range constDataTwo() const
    {
        if (d->size == 0 || isLinearised())
            return const_array_range(nullptr, 0);
        return const_array_range(d->data(), d->last + 1);
    }

    bool endsWith(const T &val) const { return !isEmpty() && last() == val; }
    QCircularBuffer<T> &fill(const T &val, int number = -1);
    T &first() { Q_ASSERT(!isEmpty()); d.detach(); return d->data()[ d->first ]; }
    const T &first() const { Q_ASSERT(!isEmpty()); return d->data()[ d->first ]; }
    int freeSize() const { return sizeAvailable(); }

    static QCircularBuffer<T> fromList(const QList<T> &list)
    { return QCircularBuffer(list.begin(), list.end()); }
    static QCircularBuffer<T> fromVector(const QVector<T> &vector)
    { return QCircularBuffer(vector.begin(), vector.end()); }

    int indexOf(const T &val, int from = 0) const;
    void insert(int i, const T &val) { insert(i, 1, val); }
    void insert(int i, int number, const T &val);
    bool isEmpty() const { return d->size == 0; }
    bool isFull() const { return d->size == d->capacity; }
    bool isLinearised() const { return (d->last >= d->first); }
    T &last() { Q_ASSERT(!isEmpty()); return d->data()[ d->last ]; }
    const T &last() const { Q_ASSERT(!isEmpty()); return d->data()[ d->last ]; }
    int lastIndexOf(const T &val, int from = -1) const;
    void linearise()
    {
        if (!isLinearised()) {
            QCircularBuffer linearized(this->cbegin(), this->cend());
            swap(linearized);
        }
    }

    void prepend(const T &val);
    void remove(int i) { remove(i, 1); }
    void remove(int i, int number);

    void replace(int i, const T &val)
    {
        Q_ASSERT_X(i >= 0 && i < d->size, "QCircularBuffer<T>::replace", "index out of range");
        const T copy(val);
        (*this)[ i ] = copy;
    }

    void reserve(int amount) { setCapacity(amount); }
    void resize(int newSize);
    void setCapacity(int amount);
    int size() const { return d->size; }
    Q_DECL_CONSTEXPR int max_size() const { return std::numeric_limits<size_type>::max(); }
    int sizeAvailable() const { return d->capacity - d->size; }
    void squeeze() { setCapacity(size()); }
    bool startsWith(const T &val) const { return !isEmpty() && first() == val; }

    QList<T> toList() const;
    QVector<T> toVector() const;

    T value(int i) const
    {
        if (i < 0 || i >= d->size)
            return T();
        return at(i);
    }

    T value(int i, const T &defaultValue) const
    {
        if (i < 0 || i >= d->size)
            return defaultValue;
        return at(i);
    }

    QCircularBuffer<T> &operator+=(const T &other) { append(other); return *this; }
    QCircularBuffer<T> &operator+=(const QCircularBuffer<T> &other);
    QCircularBuffer<T> &operator+=(const QVector<T> &other);
    QCircularBuffer<T> &operator+=(const QList<T> &other);

    QCircularBuffer<T> &operator<<(const T &other) { append(other); return *this; }
    QCircularBuffer<T> &operator<<(const QCircularBuffer<T> &other) { *this += other; return *this; }
    QCircularBuffer<T> &operator<<(const QVector<T> &other) { *this += other; return *this; }
    QCircularBuffer<T> &operator<<(const QList<T> &other) { *this += other; return *this; }

    inline bool isSharedWith(const QCircularBuffer &other) const { return d == other.d; }

private:
    QExplicitlySharedDataPointer<Data> d;
};

template <typename T>
QCircularBuffer<T> operator+(const QCircularBuffer<T> &lhs, const QCircularBuffer<T> &rhs);

template <typename T>
inline bool operator==(const QCircularBuffer<T> &lhs, const QCircularBuffer<T> &rhs)
{ return lhs.isSharedWith(rhs) || (lhs.size() == rhs.size() && std::equal(lhs.begin(), lhs.end(), rhs.begin())); }

template <typename T>
inline bool operator!=(const QCircularBuffer<T> &lhs, const QCircularBuffer<T> &rhs)
{ return !operator==(lhs, rhs); }

template <typename T>
inline void swap(QCircularBuffer<T> &lhs, QCircularBuffer<T> &rhs)
{ lhs.swap(rhs); }

template <typename T>
inline bool operator< (const QCircularBuffer<T> &lhs, const QCircularBuffer<T> &rhs)
{ return std::lexicographical_compare(lhs.begin(), lhs.end(), rhs.begin(), rhs.end()); }

template <typename T>
inline bool operator> (const QCircularBuffer<T> &lhs, const QCircularBuffer<T> &rhs)
{ return operator<(rhs, lhs); }

template <typename T>
inline bool operator>=(const QCircularBuffer<T> &lhs, const QCircularBuffer<T> &rhs)
{ return !operator<(lhs, rhs); }

template <typename T>
inline bool operator<=(const QCircularBuffer<T> &lhs, const QCircularBuffer<T> &rhs)
{ return !operator>(lhs, rhs); }

// out-of-line function implementations:

#ifndef Q_QDOC

template <typename T>
QCircularBuffer<T>::QCircularBuffer(int amount)
    : d(new Data())
{
    // Allocate some raw memory
    d->setData(d->allocate(amount));
    d->capacity = amount;

    // Initialize memory block to zero
    memset(d->data(), 0, amount * sizeof(T));
}

template <typename T>
QCircularBuffer<T>::QCircularBuffer(int amount, const T &val)
    : d(new Data())
{
    // Allocate some raw memory
    d->setData(d->allocate(amount));
    d->capacity = amount;

    // Initialize the objects. In this case we always use the placement new operator
    T *b = d->data();
    T *i = b + d->capacity;
    while (i != b)
        new (--i) T(val);
    d->first = 0;
    d->last = d->capacity - 1;
    d->size = d->capacity;
}

template <typename T>
QCircularBuffer<T>::QCircularBuffer(int amount, int initialSize, const T &val)
    : d(new Data())
{
    Q_ASSERT_X(amount >= initialSize, "QCircularBuffer<T>::QCircularBuffer(int capacity, int size, const T &value)", "size is greater than capacity");

    // Allocate some raw memory
    d->setData(d->allocate(amount));
    d->capacity = amount;

    // Initialize the objects that need to be set to the specified value.
    // In this case we always use the placement new operator
    T *b = d->data();
    T *i = b + initialSize;
    while (i != b)
        new (--i) T(val);

    // Initialize the remaining memory to zero
    memset(d->data() + initialSize, 0, (amount - initialSize) * sizeof(T));

    d->first = 0;
    d->last = initialSize - 1;
    d->size = initialSize;
}

template <typename T>
void QCircularBuffer<T>::append(const T &val)
{
    // If we have no capacity we do nothing
    if (!d->capacity)
        return;
    d.detach();
    if (d->size == d->capacity) {
        // Buffer is full. Overwrite earliest item and rotate
        d->data()[ d->first ] = val;
        d->first = d->wraparound(++d->first);
        d->last = d->wraparound(++d->last);
    } else if (d->size != 0) {
        // Buffer is partially full. Append data to end of array using appropriate method
        int index = d->wraparound(d->first + d->size);
        if (QTypeInfo<T>::isComplex)
            new (d->data() + index) T(val);
        else
            d->data()[ index ] = val;
        ++d->size;
        ++d->last;
    } else {
        // Buffer is empty. Append data to end of array using appropriate method
        d->size = 1;
        d->first = d->last = 0;
        if (QTypeInfo<T>::isComplex)
            new (d->data() + d->first) T(val);
        else
            d->data()[ d->first ] = val;
    }
}

template <typename T>
bool QCircularBuffer<T>::contains(const T &val) const
{
    if (isLinearised()) {
        T *b = d->data() + d->first;
        T *i = b + d->size;
        while (i != b)
            if (*--i == val)
                return true;
        return false;
    } else {
        // Check the array from m_first to the end
        T *b = d->data() + d->first;
        T *i = d->data() + d->capacity;
        while (i != b)
            if (*--i == val)
                return true;

        // Check array from 0 to m_end
        b = d->data();
        i = d->data() + d->last + 1;
        while (i != b)
            if (*--i == val)
                return true;

        return false;
    }
}

template <typename T>
int QCircularBuffer<T>::count(const T &val) const
{
    int c = 0;
    if (isLinearised()) {
        T *b = d->data() + d->first;
        T *i = b + d->size;
        while (i != b)
            if (*--i == val)
                ++c;
    } else {
        // Check the array from m_first to the end
        T *b = d->data() + d->first;
        T *i = d->data() + d->capacity;
        while (i != b)
            if (*--i == val)
                ++c;

        // Check array from 0 to m_end
        b = d->data();
        i = d->data() + d->last + 1;
        while (i != b)
            if (*--i == val)
                ++c;
    }
    return c;
}

template <typename T>
QCircularBuffer<T> &QCircularBuffer<T>::fill(const T &val, int number)
{
    Q_ASSERT_X(d->capacity >= number, "QCircularBuffer<T>::fill", "size is greater than capacity");
    const T copy(val);
    d.detach();
    int oldSize = d->size;
    d->size = (number < 0 ? d->size : number);
    d->first = (number == 0 ? -1 : 0);
    d->last = d->size - 1;

    // Copy item into array size times
    if (d->size) {
        T *b = d->data();
        T *i = d->data() + d->size;
        while (i != b)
            *--i = copy;
    }

    if (d->size < oldSize) {
        // Cleanup old items beyond end of new array
        T *b = d->data() + d->size;
        T *i = d->data() + oldSize;
        while (i-- != b) {
            i->~T();
            //TODO: Optimize to a single memset call
            memset(i, 0, sizeof(T));
        }
    }

    return *this;
}

template <typename T>
int QCircularBuffer<T>::indexOf(const T &val, int from) const
{
    Q_ASSERT_X(from < d->size, "QCircularBuffer<T>::indexOf", "from is greater than last valid index");
    if (from < 0)
        from = qMax(from + d->size, 0);
    else if (from >= d->size)
        from = d->size - 1;
    for (int i = from; i < d->size; ++i)
        if (at(i) == val)
            return i;
    return -1;
}

template <typename T>
void QCircularBuffer<T>::insert(int i, int number, const T &val)
{
    Q_ASSERT_X(i >= 0 && i <= d->size, "QCircularBuffer<T>::insert", "index out of range");
    d.detach();
    int freeCapacity = d->capacity - d->size;

    // Calculate number of elements that will actually be inserted. This
    // depends upon where the insertion has been requested and any spare
    // capacity left in the buffer. This is because elements at higher
    // indices will be pushed to the right and will potentially wrap around
    // to overwrite earlier elements.
    int numToInsert = qMin(number, i + freeCapacity);

    // Calculate the number of elements at the beginning of the buffer that
    // will be overwritten
    int numToOverwrite = qMin(i, qMax(0, number - freeCapacity));

    // Decide which way to shift to minimize the amount of copying required.
    if (i < d->size / 2) {
        // Inserting in lower half of buffer so we shift earlier items down

        // Shift data at the bottom end down. This may only be a subset if some
        // of the early data is to be overwritten.
        if (QTypeInfo<T>::isStatic) {
            int b = d->first + numToOverwrite;
            int e = d->first + i - 1;
            for (int j = b; j <= e; ++j) {
                int srcIndex = j % d->capacity;
                int dstIndex = (j - numToInsert + d->capacity) % d->capacity;
                T *src = d->data() + srcIndex;
                T *dst = d->data() + dstIndex;

                new (dst) T(*src);
            }
        } else {
            // We have a movable type so a simple memcopy (or maybe two or
            // three) will suffice to shift the data at the bottom end
            int numToMove = i - numToOverwrite;
            if (numToMove > 0) {
                int srcBegin = (d->first + numToOverwrite) % d->capacity;
                int srcEnd = (d->first + i - 1) % d->capacity;
                int dstBegin = (srcBegin - numToInsert + d->capacity) % d->capacity;
                int dstEnd = (srcEnd - numToInsert + d->capacity) % d->capacity;

                // Do we have any wrap-around problems to deal with?
                bool srcRegionWraps = (srcEnd < srcBegin);
                bool dstRegionWraps = (dstEnd < dstBegin);
                if (!srcRegionWraps && dstRegionWraps) {
                    // Destination region wraps so do the move in two steps
                    int wrapCount = abs(srcBegin - numToInsert);
                    memmove(d->data() + d->capacity - wrapCount, d->data() + srcBegin, wrapCount * sizeof(T));
                    memmove(d->data(), d->data() + srcBegin + wrapCount, (numToMove - wrapCount) * sizeof(T));
                } else if (srcRegionWraps && !dstRegionWraps) {
                    // Source region wraps so do the move in two steps
                    int wrapCount = d->capacity - srcBegin;
                    memmove(d->data() + dstBegin, d->data() + d->capacity - wrapCount, wrapCount * sizeof(T));
                    memmove(d->data() + dstBegin + numToInsert, d->data(), (numToMove - wrapCount) * sizeof(T));
                } else if (srcRegionWraps && dstRegionWraps) {
                    // Source and destination regions wrap so we have to do this in three steps
                    int srcWrapCount = d->capacity - srcBegin;
                    memmove(d->data() + dstBegin, d->data() + d->capacity - srcWrapCount, srcWrapCount * sizeof(T));
                    memmove(d->data() + d->capacity - numToInsert, d->data(), numToInsert * sizeof(T));
                    memmove(d->data(), d->data() + numToInsert, (numToMove - srcWrapCount - numToInsert) * sizeof(T));
                } else {
                    // No wrap around - do a single memmove
                    memmove(d->data() + dstBegin, d->data() + srcBegin, numToMove * sizeof(T));
                }
            }
        }

        // Insert the new items
        int e = d->first + i;
        int b = e - numToInsert;
        for (int j = b; j < e; ++j) {
            T *p = d->data() + ((j + d->capacity) % d->capacity);
            new (p) T(val);
        }

        // Adjust the first, last and size indices as needed.
        // NB. The last index never changes in this regime.
        d->size += qMin(number, freeCapacity);
        d->first = (d->first - (numToInsert - numToOverwrite) + d->capacity) % d->capacity;
    } else {
        // Inserting in upper half of buffer so we shift later items up

        // Shift data at the top end up which may or may not overwrite some
        // of the earliest data
        if (QTypeInfo<T>::isStatic) {
            int b = d->first + d->size - 1;
            int e = d->first + i;
            for (int j = b; j >= e; j--) {
                int srcIndex = j % d->capacity;
                int dstIndex = (j + numToInsert) % d->capacity;
                T *src = d->data() + srcIndex;
                T *dst = d->data() + dstIndex;

                new (dst) T(*src);
            }
        } else {
            // We have a movable type so a simple memcopy (or maybe two or
            // three) will suffice to shift the data at the top end
            int numToMove = d->size - i;
            if (numToMove > 0) {
                int srcBegin = (d->first + i) % d->capacity;
                int srcEnd = d->last;
                int dstBegin = (srcBegin + numToInsert) % d->capacity;
                int dstEnd = (srcEnd + numToInsert) % d->capacity;

                // Do we have any wrap-around problems to deal with?
                bool srcRegionWraps = (srcEnd < srcBegin);
                bool dstRegionWraps = (dstEnd < dstBegin);
                if (!srcRegionWraps && dstRegionWraps) {
                    // Destination region wraps so do the move in two steps
                    int wrapCount = srcEnd + numToInsert - d->capacity + 1;
                    memmove(d->data(), d->data() + srcEnd - wrapCount + 1, wrapCount * sizeof(T));
                    memmove(d->data() + dstBegin, d->data() + srcBegin, (numToMove - wrapCount) * sizeof(T));
                } else if (srcRegionWraps && !dstRegionWraps) {
                    // Source region wraps so do the move in two steps
                    int wrapCount = d->last + 1;
                    memmove(d->data() + numToInsert, d->data(), wrapCount * sizeof(T));
                    memmove(d->data() + dstBegin, d->data() + srcBegin, (numToMove - wrapCount) * sizeof(T));
                } else if (srcRegionWraps && dstRegionWraps) {
                    // Source and destination regions wrap so we have to do this in three steps
                    int srcWrapCount = d->last + 1;
                    memmove(d->data() + numToInsert, d->data(), srcWrapCount * sizeof(T));
                    memmove(d->data(), d->data() + d->capacity - numToInsert, numToInsert * sizeof(T));
                    memmove(d->data() + dstBegin, d->data() + srcBegin, (numToMove - srcWrapCount - numToInsert) * sizeof(T));
                } else {
                    // No wrap around - do a single memmove
                    memmove(d->data() + dstBegin, d->data() + srcBegin, numToMove * sizeof(T));
                }
            }
        }

        // Insert the new items
        for (int j = d->first + i; j < d->first + i + numToInsert; ++j) {
            T *p = d->data() + (j % d->capacity);
            new (p) T(val);
        }

        // Adjust the first, last and size indices as needed
        d->size += qMin(number, freeCapacity);
        d->first = (d->first + numToOverwrite) % d->capacity;
        d->last = (d->last + numToInsert) % d->capacity;
    }
}

template <typename T>
int QCircularBuffer<T>::lastIndexOf(const T &val, int from) const
{
    if (from < 0)
        from = qMax(from + d->size, 0);
    else if (from >= d->size)
        from = d->size - 1;
    for (int i = from; i >= 0; --i)
        if (at(i) == val)
            return i;
    return -1;
}

template <typename T>
void QCircularBuffer<T>::prepend(const T &val)
{
    // If we have no capacity we do nothing
    if (!d->capacity)
        return;

    d.detach();
    if (d->size == d->capacity) {
        // Buffer is full. Overwrite last item and rotate
        d->data()[ d->last ] = val;
        d->first = (--d->first + d->capacity) % d->capacity;
        d->last = (--d->last + d->capacity) % d->capacity;
    } else if (d->size != 0) {
        // Buffer is partially full. Prepend data to start of array using appropriate method
        d->first = (--d->first + d->capacity) % d->capacity;
        ++d->size;
        if (QTypeInfo<T>::isComplex)
            new (d->data() + d->first) T(val);
        else
            d->data()[ d->first ] = val;
    } else {
        // Buffer is empty. Prepend data to start of array using appropriate method
        d->size = 1;
        d->first = d->last = d->capacity - 1;
        if (QTypeInfo<T>::isComplex)
            new (d->data() + d->first) T(val);
        else
            d->data()[ d->first ] = val;
    }
}

template <typename T>
void QCircularBuffer<T>::remove(int i, int number)
{
    Q_ASSERT_X(i >= 0 && number > 0 && i + number <= d->size, "QCircularBuffer<T>::remove", "index out of range");
    d.detach();

    // HACK (it actually makes sense, but requires some more thinking)
    if ( i == 0 && !QTypeInfo<T>::isComplex ) {
        d->first = d->wraparound( d->first + number );
        d->size -= number;
        return;
    }

    // Calculate the number of items that need to be moved downward
    int numToMoveDown = d->size - number - i;
    int numToMoveUp = i;

    if (numToMoveDown < numToMoveUp) {
        // Move higher items down
        int numToMove = numToMoveDown;

        if (QTypeInfo<T>::isComplex) {
            // Copy items down from higher positions
            int b = d->first + i;
            int e = b + numToMove;
            for (int j = b; j < e ; ++j) {
                T *src = d->data() + ((j + number) % d->capacity);
                T *dst = d->data() + (j % d->capacity);
                new (dst) T(*src);
            }

            // Clean up items at end of buffer
            for (int j = d->last; j > d->last - number; --j) {
                T *p = d->data() + ((j + d->capacity) % d->capacity);
                p->~T();
                //TODO: Optimize to a single memset call
                memset(p, 0, sizeof(T));
            }
        } else {
            if (isLinearised()) {
                // With a linearised buffer we can do a simple move and removal of items
                memmove(d->data() + d->last - numToMove - number + 1, d->data() + d->last - numToMove + 1, numToMove * sizeof(T));
                memset(d->data() + d->last - number + 1, 0, number * sizeof(T));
            } else {
                // With a non-linearised buffer we need to be careful of wrapping issues
                int srcBegin = (d->last - numToMove + 1 + d->capacity) % d->capacity;
                int srcEnd = d->last;
                int dstBegin = (d->first + i) % d->capacity;
                int dstEnd = (dstBegin + numToMove - 1) % d->capacity;

                bool srcRegionWraps = (srcEnd < srcBegin);
                bool dstRegionWraps = (dstEnd < dstBegin);
                if (srcRegionWraps && !dstRegionWraps) {
                    // Source region wraps so do the move in two steps
                    int wrapCount = d->capacity - srcBegin;
                    memmove(d->data() + dstBegin, d->data() + srcBegin, wrapCount * sizeof(T));
                    memmove(d->data() + dstBegin + wrapCount, d->data(), (numToMove - wrapCount) * sizeof(T));
                } else if (!srcRegionWraps && dstRegionWraps) {
                    // Destination region wraps so do the move in two steps
                    int wrapCount = number - srcBegin;
                    memmove(d->data() + d->capacity - wrapCount, d->data() + srcBegin, wrapCount * sizeof(T));
                    memmove(d->data(), d->data() + srcBegin + wrapCount, (numToMove - wrapCount) * sizeof(T));
                } else if (srcRegionWraps && dstRegionWraps) {
                    // Source and destination regions wrap so we have to do this in three steps
                    int srcWrapCount = d->capacity - srcBegin;
                    memmove(d->data() + dstBegin, d->data() + srcBegin, srcWrapCount * sizeof(T));
                    memmove(d->data() + dstBegin + srcWrapCount, d->data(), number * sizeof(T));
                    memmove(d->data(), d->data() + number, (numToMove - srcWrapCount - number) * sizeof(T));
                } else {
                    // No wrap around, so we can do this in one hit
                    memmove(d->data() + dstBegin, d->data() + srcBegin, numToMove * sizeof(T));
                }

                // We potentially have a disjoint region that needs zeroing
                int zeroStart = (d->last - number + d->capacity + 1) % d->capacity;
                int zeroEnd = d->last;
                if (zeroEnd < zeroStart) {
                    // Region to be zeroed wraps. Do it in two steps.
                    memset(d->data(), 0, (d->last + 1) * sizeof(T));
                    memset(d->data() + zeroStart, 0, (number - d->last - 1) * sizeof(T));
                } else {
                    // Region to be zeroed is contiguous
                    memset(d->data() + zeroStart, 0, number * sizeof(T));
                }
            }
        }

        // Adjust the indices
        d->size -= number;
        d->last = (d->last - number + d->capacity) % d->capacity;
    } else {
        // Move lower items up
        int numToMove = numToMoveUp;

        if (QTypeInfo<T>::isComplex) {
            // Copy items up from lower positions
            int b = d->first + i - 1;
            int e = b - numToMove;
            for (int j = b; j > e ; --j) {
                T *src = d->data() + ((j + d->capacity) % d->capacity);
                T *dst = d->data() + ((j + d->capacity + number) % d->capacity);
                new (dst) T(*src);
            }

            // Clean up items at start of buffer
            for (int j = d->first; j < d->first + number; ++j) {
                T *p = d->data() + (j % d->capacity);
                p->~T();
                //TODO: Optimize to a single memset call
                memset(p, 0, sizeof(T));
            }
        } else {
            if (isLinearised()) {
                // With a linearised buffer we can do a simple move and removal of items
                memmove(d->data() + d->first + number, d->data() + d->first, numToMove * sizeof(T));
                memset(d->data() + d->first, 0, number * sizeof(T));
            } else {
                // With a non-linearised buffer we need to be careful of wrapping issues
                int srcBegin = d->first;
                int srcEnd = (srcBegin + numToMove - 1) % d->capacity;
                int dstBegin = (srcBegin + number) % d->capacity;
                int dstEnd = (dstBegin + numToMove - 1) % d->capacity;

                bool srcRegionWraps = (srcEnd < srcBegin);
                bool dstRegionWraps = (dstEnd < dstBegin);
                if (srcRegionWraps && !dstRegionWraps) {
                    // Source region wraps so do the move in two steps
                    int wrapCount = srcEnd + 1;
                    memmove(d->data() + dstEnd - wrapCount + 1, d->data(), wrapCount * sizeof(T));
                    memmove(d->data() + dstBegin, d->data() + srcBegin, (numToMove - wrapCount) * sizeof(T));
                } else if (!srcRegionWraps && dstRegionWraps) {
                    // Destination region wraps so do the move in two steps
                    int wrapCount = dstEnd + 1;
                    memmove(d->data(), d->data() + srcEnd - wrapCount + 1, wrapCount * sizeof(T));
                    memmove(d->data() + dstBegin, d->data() + srcBegin, (numToMove - wrapCount) * sizeof(T));
                } else if (srcRegionWraps && dstRegionWraps) {
                    // Source and destination regions wrap so we have to do this in three steps
                    int srcWrapCount = srcEnd + 1;
                    memmove(d->data() + dstEnd - srcWrapCount + 1, d->data(), srcWrapCount * sizeof(T));
                    memmove(d->data(), d->data() + d->capacity - number, number * sizeof(T));
                    memmove(d->data() + dstBegin, d->data() + srcBegin, (numToMove - srcWrapCount - number) * sizeof(T));
                } else {
                    // No wrap around, so we can do this in one hit
                    memmove(d->data() + dstBegin, d->data() + srcBegin, numToMove * sizeof(T));
                }

                // We potentially have a disjoint region that needs zeroing
                int zeroStart = d->first;
                int zeroEnd = (zeroStart + number - 1) % d->capacity;
                if (zeroEnd < zeroStart) {
                    // Region to be zeroed wraps. Do it in two steps.
                    memset(d->data() + zeroStart, 0, (d->capacity - d->first) * sizeof(T));
                    memset(d->data(), 0, (number - d->capacity + d->first) * sizeof(T));
                } else {
                    // Region to be zeroed is contiguous
                    memset(d->data() + zeroStart, 0, number * sizeof(T));
                }
            }
        }

        // Adjust the indices
        d->size -= number;
        d->first = (d->first + number) % d->capacity;
    }
}

template <typename T>
void QCircularBuffer<T>::setCapacity(int amount)
{
    if (amount == d->capacity)
        return;

    d.detach();
    // Allocate some new raw memory
    T *newData = d->allocate(amount);

    // How many items can we copy across?
    int newSize = qMin(d->size, amount);

    if (QTypeInfo<T>::isComplex) {
        // Copy across the elements from the original array
        for (int i = 0; i < newSize; ++i) {
            T *src = d->data() + ((d->first + i) % d->capacity);
            T *dst = newData + i;
            new (dst) T(*src);
        }

        // Destroy the original items.
        // The type is complex so we manually call the destructor for each item
        // since we used the placement new operator to instantiate them
        T *b = d->data();
        T *i = b + d->capacity;
        while (i-- != b)
             i->~T();
    } else {
        // Copy across the elements from the original array. The source region
        // potentially wraps so we may have to do this in one or two steps
        if (isLinearised()) {
            memmove(newData, d->data() + d->first, newSize * sizeof(T));
        } else {
            int step1Size = qMin(newSize, d->capacity - d->first);
            memmove(newData, d->data() + d->first, step1Size * sizeof(T));
            int step2Size = qMax(0, qMin(newSize - d->capacity + d->first, d->last + 1));
            memmove(newData + step1Size, d->data(), step2Size * sizeof(T));
        }
    }

    // Initialize any memory outside of the valid buffer (ie the unused items)
    memset(newData + newSize, 0, (amount - newSize) * sizeof(T));

    // Release the raw memory for the old array
    d->deallocate(d->data());

    // Assign the new memory to be our buffer data and fix indices
    d->setData(newData);
    d->capacity = amount;
    d->first = 0;
    d->size = newSize;
    d->last = d->size - 1;
}

template <typename T>
void QCircularBuffer<T>::resize(int newSize)
{
    Q_ASSERT_X(newSize >= 0 && newSize <= d->capacity, "QCircularBuffer<T>::resize", "size out of range");
    d.detach();
    if (newSize < d->size) {
        remove(newSize, d->size - newSize);
    } else if (newSize > d->size) {
        const T t = T();
        insert(d->size, newSize - d->size, t);
    }
}

template <typename T>
QCircularBuffer<T> &QCircularBuffer<T>::operator+=(const QCircularBuffer<T> &other)
{
    d.detach();
    // How many items do we need to copy? No point in ever copying across a number
    // greater than capacity
    int numToCopy = qMin(other.size(), d->capacity);
    int offset = other.size() - numToCopy;
    for (int i = 0; i < numToCopy; ++i)
        append(other.at(offset + i));
    return *this;
}

template <typename T>
QCircularBuffer<T> &QCircularBuffer<T>::operator+=(const QVector<T> &other)
{
    d.detach();
    // How many items do we need to copy? No point in ever copying across a number
    // greater than capacity
    int numToCopy = qMin(other.size(), d->capacity);
    int offset = other.size() - numToCopy;
    for (int i = 0; i < numToCopy; ++i)
        append(other.at(offset + i));
    return *this;
}

template <typename T>
QCircularBuffer<T> &QCircularBuffer<T>::operator+=(const QList<T> &other)
{
    d.detach();
    // How many items do we need to copy? No point in ever copying across a number
    // greater than capacity
    int numToCopy = qMin(other.size(), d->capacity);
    int offset = other.size() - numToCopy;
    for (int i = 0; i < numToCopy; ++i)
        append(other.at(offset + i));
    return *this;
}

template <typename T>
QList<T> QCircularBuffer<T>::toList() const
{
    QList<T> list;
    list.reserve(size());
    for (int i = 0; i < size(); ++i)
        list.append(at(i));
    return list;
}

template <typename T>
QVector<T> QCircularBuffer<T>::toVector() const
{
    QVector<T> vector;
    vector.reserve(size());
    for (int i = 0; i < size(); ++i)
        vector.append(at(i));
    return vector;
}

template <typename T>
QCircularBuffer<T> operator+(const QCircularBuffer<T> &lhs, const QCircularBuffer<T> &rhs)
{
    QCircularBuffer<T> circ(lhs.size() + rhs.size());
    for (int i = 0; i < lhs.size(); ++i)
        circ.append(lhs.at(i));
    for (int i = 0; i < rhs.size(); ++i)
        circ.append(rhs.at(i));
    return circ;
}

#endif // Q_QDOC

Q_DECLARE_SEQUENTIAL_ITERATOR(CircularBuffer)
Q_DECLARE_MUTABLE_SEQUENTIAL_ITERATOR(CircularBuffer)

} //Qt3D

QT_END_NAMESPACE

#endif // QT3DCORE_QCIRCULARBUFFER_H
