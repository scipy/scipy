/****************************************************************************
**
** Copyright (C) 2008-2012 NVIDIA Corporation.
** Copyright (C) 2019 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of Qt Quick 3D.
**
** $QT_BEGIN_LICENSE:GPL$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 3 or (at your option) any later version
** approved by the KDE Free Qt Foundation. The licenses are as published by
** the Free Software Foundation and appearing in the file LICENSE.GPL3
** included in the packaging of this file. Please review the following
** information to ensure the GNU General Public License requirements will
** be met: https://www.gnu.org/licenses/gpl-3.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/


#ifndef QSSGPERFRAMEALLOCATOR_H
#define QSSGPERFRAMEALLOCATOR_H

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API.  It exists purely as an
// implementation detail.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.
//

#include <QtQuick3DRuntimeRender/private/qtquick3druntimerenderglobal_p.h>

QT_BEGIN_NAMESPACE

class QSSGPerFrameAllocator
{
    struct FastAllocator
    {
        struct Slab;

        enum : size_t {
            ChunkSize = 8192*2,
            Alignment = 4,
            SlabSize = ChunkSize - sizeof(Slab *),
            MaxAlloc = ChunkSize/2 // don't go all the way up to SlabSize, or we'd almost always get a big hole
        };
        struct Slab {
            Slab() = default;
            Slab(Slab *previous)
            {
                previous->next = this;
            }
            Slab *next = nullptr;
            quint8 data[SlabSize];
        };
        Q_STATIC_ASSERT(sizeof(Slab) == ChunkSize);

        Slab *first = nullptr;
        Slab *current = nullptr;
        size_t offset = 0;

        FastAllocator()
        {
            first = current = new Slab;
        }

        ~FastAllocator()
        {
            Slab *s = first;
            while (s) {
                Slab *n = s->next;
                delete s;
                s = n;
            }
        }
        void *allocate(size_t size)
        {
            size = (size + Alignment - 1) & ~(Alignment - 1);
            Q_ASSERT(size <= SlabSize);
            Q_ASSERT(!(offset % Alignment));

            size_t amountLeftInSlab = SlabSize - offset;
            if (size > amountLeftInSlab) {
                if (current->next)
                    current = current->next;
                else
                    current = new Slab(current);
                offset = 0;
            }

            quint8 *data = current->data + offset;
            offset += size;
            return data;
        }

        // only reset, so we can re-use the memory
        void reset() { current = first; offset = 0; }
    };

    struct LargeAllocator
    {
        struct Slab {
            Slab *next = nullptr;
        };
        Slab *current = nullptr;

        LargeAllocator() {}

        // Automatically deallocates everything that hasn't already been deallocated.
        ~LargeAllocator() { deallocateAll(); }

        void deallocateAll()
        {
            while (current) {
                Slab *n = current->next;
                ::free(current);
                current = n;
            }
            current = nullptr;
        }

        void *allocate(size_t size)
        {
            quint8 *mem = reinterpret_cast<quint8 *>(::malloc(sizeof(Slab) + size));
            Slab *s = reinterpret_cast<Slab *>(mem);
            s->next = current;
            current = s;
            return mem + sizeof(Slab);
        }
    };

    FastAllocator m_fastAllocator;
    LargeAllocator m_largeAllocator;

public:
    QSSGPerFrameAllocator() {}

    inline void *allocate(size_t size)
    {
        if (size < FastAllocator::MaxAlloc)
            return m_fastAllocator.allocate(size);

        return m_largeAllocator.allocate(size);
    }

    void reset()
    {
        m_fastAllocator.reset();
        m_largeAllocator.deallocateAll();
    }
};

QT_END_NAMESPACE

#endif // QSSGPERFRAMEALLOCATOR_H
