/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt Data Visualization module of the Qt Toolkit.
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

//
//  W A R N I N G
//  -------------
//
// This file is not part of the QtDataVisualization API.  It exists purely as an
// implementation detail.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.

#ifndef SCATTERSERIESRENDERCACHE_P_H
#define SCATTERSERIESRENDERCACHE_P_H

#include "datavisualizationglobal_p.h"
#include "seriesrendercache_p.h"
#include "qscatter3dseries_p.h"
#include "scatterrenderitem_p.h"

QT_BEGIN_NAMESPACE_DATAVISUALIZATION

class ScatterObjectBufferHelper;
class ScatterPointBufferHelper;

class ScatterSeriesRenderCache : public SeriesRenderCache
{
public:
    ScatterSeriesRenderCache(QAbstract3DSeries *series, Abstract3DRenderer *renderer);
    virtual ~ScatterSeriesRenderCache();

    void cleanup(TextureHelper *texHelper);

    inline ScatterRenderItemArray &renderArray() { return m_renderArray; }
    inline QScatter3DSeries *series() const { return static_cast<QScatter3DSeries *>(m_series); }
    inline void setItemSize(float size) { m_itemSize = size; }
    inline float itemSize() const { return m_itemSize; }
    inline void setSelectionIndexOffset(int offset) { m_selectionIndexOffset = offset; }
    inline int selectionIndexOffset() const { return m_selectionIndexOffset; }
    inline void setStaticBufferDirty(bool state) { m_staticBufferDirty = state; }
    inline bool staticBufferDirty() const { return m_staticBufferDirty; }
    inline int oldArraySize() const { return m_oldRenderArraySize; }
    inline void setOldArraySize(int size) { m_oldRenderArraySize = size; }
    inline const QString &oldMeshFileName() const { return m_oldMeshFileName; }
    inline void setOldMeshFileName(const QString &meshFileName) { m_oldMeshFileName = meshFileName; }
    inline void setBufferObject(ScatterObjectBufferHelper *object) { m_scatterBufferObj = object; }
    inline ScatterObjectBufferHelper *bufferObject() const { return m_scatterBufferObj; }
    inline void setBufferPoints(ScatterPointBufferHelper *object) { m_scatterBufferPoints = object; }
    inline ScatterPointBufferHelper *bufferPoints() const { return m_scatterBufferPoints; }
    inline QVector<int> &updateIndices() { return m_updateIndices; }
    inline QVector<int> &bufferIndices() { return m_bufferIndices; }
    inline void setVisibilityChanged(bool changed) { m_visibilityChanged = changed; }
    inline bool visibilityChanged() const { return m_visibilityChanged; }

protected:
    ScatterRenderItemArray m_renderArray;
    float m_itemSize;
    int m_selectionIndexOffset; // Temporarily cached value for selection color calculations
    bool m_staticBufferDirty;
    int m_oldRenderArraySize; // Used to detect if full buffer change needed
    QString m_oldMeshFileName; // Used to detect if full buffer change needed
    ScatterObjectBufferHelper *m_scatterBufferObj;
    ScatterPointBufferHelper *m_scatterBufferPoints;
    QVector<int> m_updateIndices; // Used as temporary cache during item updates
    QVector<int> m_bufferIndices; // Cache for mapping renderarray to mesh buffer
    bool m_visibilityChanged; // Used to detect if full buffer change needed
};

QT_END_NAMESPACE_DATAVISUALIZATION

#endif
