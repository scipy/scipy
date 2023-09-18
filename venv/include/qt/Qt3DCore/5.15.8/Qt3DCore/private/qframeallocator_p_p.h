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

#ifndef QT3DCORE_QFRAMEALLOCATOR_P_P_H
#define QT3DCORE_QFRAMEALLOCATOR_P_P_H

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

QT_BEGIN_NAMESPACE

namespace Qt3DCore {

class QFrameAllocator;

struct Q_AUTOTEST_EXPORT QFrameChunk
{
    void init(uint blockSize, uchar blocks);
    void *allocate(uint blockSize);

    void deallocate(void *p, uint blockSize);
    bool contains(void *p, uint blockSize);
    void clear(uint blockSize, uchar blocks);
    void release();

    inline bool isEmpty() const { return m_blocksAvailable == m_maxBlocksAvailable; }

    uchar *m_data;
    uchar  m_firstAvailableBlock;
    uchar  m_blocksAvailable;
    uchar  m_maxBlocksAvailable;
};

class Q_AUTOTEST_EXPORT QFixedFrameAllocator
{
public:
    QFixedFrameAllocator();
    ~QFixedFrameAllocator();

    void init(uint blockSize, uchar pageSize = 128);
    void *allocate();
    void  deallocate(void *ptr);
    void trim();
    void release();
    void clear();
    bool isEmpty() const;

    inline int chunkCount() const { return m_chunks.size(); }
    inline uchar pageSize() const { return m_nbrBlock; }
    inline uint blockSize() const { return m_blockSize; }

private:
    QFrameChunk &scan();

private:
    uint m_blockSize;
    uchar m_nbrBlock;
    QVector<QFrameChunk> m_chunks;
    QFrameChunk *m_lastAllocatedChunck;
    QFrameChunk *m_lastFreedChunck;
};

class QFrameAllocatorPrivate
{
public:
    QFrameAllocatorPrivate();

    inline void *allocateAtChunk(uint allocatorIndex)
    {
        return m_allocatorPool[allocatorIndex].allocate();
    }

    inline void deallocateAtChunck(void *ptr, uint allocatorIndex)
    {
        m_allocatorPool[allocatorIndex].deallocate(ptr);
    }

    inline uint allocatorIndexFromSize(uint targetSize) const
    {
        return (targetSize + m_alignment - 1) / m_alignment - 1;
    }

    uint m_maxObjectSize;
    uint m_alignment;
    QVector<QFixedFrameAllocator> m_allocatorPool;
};

} // Qt3D

QT_END_NAMESPACE

#endif // QT3DCORE_QFRAMEALLOCATOR_P_P_H
