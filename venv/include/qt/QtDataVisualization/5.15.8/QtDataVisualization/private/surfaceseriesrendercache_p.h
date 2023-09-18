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

#ifndef SURFACESERIESRENDERCACHE_P_H
#define SURFACESERIESRENDERCACHE_P_H

#include "datavisualizationglobal_p.h"
#include "seriesrendercache_p.h"
#include "qsurface3dseries_p.h"
#include "surfaceobject_p.h"
#include "selectionpointer_p.h"

#include <QtGui/QMatrix4x4>

QT_BEGIN_NAMESPACE_DATAVISUALIZATION

class Surface3DRenderer;

class SurfaceSeriesRenderCache : public SeriesRenderCache
{
public:
    SurfaceSeriesRenderCache(QAbstract3DSeries *series, Surface3DRenderer *renderer);
    virtual ~SurfaceSeriesRenderCache();

    virtual void populate(bool newSeries);
    virtual void cleanup(TextureHelper *texHelper);

    inline bool surfaceVisible() const { return m_surfaceVisible; }
    inline bool surfaceGridVisible() const { return m_surfaceGridVisible; }
    inline bool isFlatShadingEnabled() const { return m_surfaceFlatShading; }
    inline void setFlatShadingEnabled(bool enabled) { m_surfaceFlatShading = enabled; }
    inline void setFlatChangeAllowed(bool allowed) { m_flatChangeAllowed = allowed; }
    inline SurfaceObject *surfaceObject() { return m_surfaceObj; }
    inline SurfaceObject *sliceSurfaceObject() { return m_sliceSurfaceObj; }
    inline const QRect &sampleSpace() const { return m_sampleSpace; }
    inline void setSampleSpace(const QRect &sampleSpace) { m_sampleSpace = sampleSpace; }
    inline QSurface3DSeries *series() const { return static_cast<QSurface3DSeries *>(m_series); }
    inline QSurfaceDataArray &dataArray() { return m_dataArray; }
    inline QSurfaceDataArray &sliceDataArray() { return m_sliceDataArray; }
    inline bool renderable() const { return m_visible && (m_surfaceVisible ||
                                                          m_surfaceGridVisible); }
    inline void setSelectionTexture(GLuint texture) { m_selectionTexture = texture; }
    inline GLuint selectionTexture() const { return m_selectionTexture; }
    inline void setSelectionIdRange(uint start, uint end) { m_selectionIdStart = start;
                                                            m_selectionIdEnd = end; }
    inline uint selectionIdStart() const { return m_selectionIdStart; }
    inline bool isWithinIdRange(uint selection) const { return selection >= m_selectionIdStart &&
                                                        selection <= m_selectionIdEnd; }
    inline bool isFlatStatusDirty() const { return m_flatStatusDirty; }
    inline void setFlatStatusDirty(bool status) { m_flatStatusDirty = status; }
    inline void setMVPMatrix(const QMatrix4x4 &matrix) { m_MVPMatrix = matrix; }
    inline const QMatrix4x4 &MVPMatrix() { return m_MVPMatrix; }

    inline void setSliceSelectionPointer(SelectionPointer *pointer) { m_sliceSelectionPointer = pointer; }
    inline SelectionPointer *sliceSelectionPointer() const { return m_sliceSelectionPointer; }
    inline void setMainSelectionPointer(SelectionPointer *pointer) { m_mainSelectionPointer = pointer; }
    inline SelectionPointer *mainSelectionPointer() const { return m_mainSelectionPointer; }

    inline void setSlicePointerActivity(bool activity) { m_slicePointerActive = activity; }
    inline bool slicePointerActive() const { return m_slicePointerActive; }
    inline void setMainPointerActivity(bool activity) { m_mainPointerActive = activity; }
    inline bool mainPointerActive() const { return m_mainPointerActive; }
    inline void setSurfaceTexture(GLuint texture) { m_surfaceTexture = texture; }
    inline GLuint surfaceTexture() const { return m_surfaceTexture; }

protected:
    bool m_surfaceVisible;
    bool m_surfaceGridVisible;
    bool m_surfaceFlatShading;
    SurfaceObject *m_surfaceObj;
    SurfaceObject *m_sliceSurfaceObj;
    QRect m_sampleSpace;
    QSurfaceDataArray m_dataArray;
    QSurfaceDataArray m_sliceDataArray;
    GLuint m_selectionTexture;
    uint m_selectionIdStart;
    uint m_selectionIdEnd;
    bool m_flatChangeAllowed;
    bool m_flatStatusDirty;
    QMatrix4x4 m_MVPMatrix;
    SelectionPointer *m_sliceSelectionPointer;
    SelectionPointer *m_mainSelectionPointer;
    bool m_slicePointerActive;
    bool m_mainPointerActive;
    GLuint m_surfaceTexture;
};

QT_END_NAMESPACE_DATAVISUALIZATION

#endif
