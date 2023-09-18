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

#ifndef SURFACEOBJECT_P_H
#define SURFACEOBJECT_P_H

#include "datavisualizationglobal_p.h"
#include "abstractobjecthelper_p.h"
#include "qsurfacedataproxy.h"

#include <QtCore/QRect>

QT_BEGIN_NAMESPACE_DATAVISUALIZATION

class Surface3DRenderer;
class AxisRenderCache;

class SurfaceObject : public AbstractObjectHelper
{
public:
    enum SurfaceType {
        SurfaceSmooth,
        SurfaceFlat,
        Undefined
    };

    enum DataDimension {
        BothAscending = 0,
        XDescending = 1,
        ZDescending = 2,
        BothDescending = XDescending | ZDescending
    };
    Q_DECLARE_FLAGS(DataDimensions, DataDimension)

public:
    SurfaceObject(Surface3DRenderer *renderer);
    virtual ~SurfaceObject();

    void setUpData(const QSurfaceDataArray &dataArray, const QRect &space,
                   bool changeGeometry, bool polar, bool flipXZ = false);
    void setUpSmoothData(const QSurfaceDataArray &dataArray, const QRect &space,
                         bool changeGeometry, bool polar, bool flipXZ = false);
    void smoothUVs(const QSurfaceDataArray &dataArray, const QSurfaceDataArray &modelArray);
    void coarseUVs(const QSurfaceDataArray &dataArray, const QSurfaceDataArray &modelArray);
    void updateCoarseRow(const QSurfaceDataArray &dataArray, int rowIndex, bool polar);
    void updateSmoothRow(const QSurfaceDataArray &dataArray, int startRow, bool polar);
    void updateSmoothItem(const QSurfaceDataArray &dataArray, int row, int column, bool polar);
    void updateCoarseItem(const QSurfaceDataArray &dataArray, int row, int column, bool polar);
    void createSmoothIndices(int x, int y, int endX, int endY);
    void createCoarseSubSection(int x, int y, int columns, int rows);
    void createSmoothGridlineIndices(int x, int y, int endX, int endY);
    void createCoarseGridlineIndices(int x, int y, int endX, int endY);
    void uploadBuffers();
    GLuint gridElementBuf();
    GLuint uvBuf();
    GLuint gridIndexCount();
    QVector3D vertexAt(int column, int row);
    void clear();
    float minYValue() const { return m_minY; }
    float maxYValue() const { return m_maxY; }
    inline void activateSurfaceTexture(bool value) { m_returnTextureBuffer = value; }

private:
    void createCoarseIndices(GLint *indices, int &p, int row, int upperRow, int j);
    void createNormals(int &p, int row, int upperRow, int j);
    void createSmoothNormalBodyLine(int &totalIndex, int column);
    void createSmoothNormalUpperLine(int &totalIndex);
    QVector3D createSmoothNormalBodyLineItem(int x, int y);
    QVector3D createSmoothNormalUpperLineItem(int x, int y);
    QVector3D normal(const QVector3D &a, const QVector3D &b, const QVector3D &c);
    void createBuffers(const QVector<QVector3D> &vertices, const QVector<QVector2D> &uvs,
                       const QVector<QVector3D> &normals, const GLint *indices);
    void checkDirections(const QSurfaceDataArray &array);
    inline void getNormalizedVertex(const QSurfaceDataItem &data, QVector3D &vertex, bool polar,
                                    bool flipXZ);

private:
    SurfaceType m_surfaceType = Undefined;
    int m_columns = 0;
    int m_rows = 0;
    GLuint m_gridElementbuffer;
    GLuint m_gridIndexCount = 0;
    QVector<QVector3D> m_vertices;
    QVector<QVector3D> m_normals;
    // Caches are not owned
    AxisRenderCache &m_axisCacheX;
    AxisRenderCache &m_axisCacheY;
    AxisRenderCache &m_axisCacheZ;
    Surface3DRenderer *m_renderer;
    float m_minY;
    float m_maxY;
    GLuint m_uvTextureBuffer;
    bool m_returnTextureBuffer = false;
    SurfaceObject::DataDimensions m_dataDimension;
    SurfaceObject::DataDimensions m_oldDataDimension = DataDimensions(-1);
};

QT_END_NAMESPACE_DATAVISUALIZATION

#endif
