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

#ifndef SURFACE3DRENDERER_P_H
#define SURFACE3DRENDERER_P_H

#include "datavisualizationglobal_p.h"
#include "surface3dcontroller_p.h"
#include "abstract3drenderer_p.h"
#include "surfaceseriesrendercache_p.h"

QT_BEGIN_NAMESPACE_DATAVISUALIZATION

class ShaderHelper;
class Q3DScene;

class QT_DATAVISUALIZATION_EXPORT Surface3DRenderer : public Abstract3DRenderer
{
    Q_OBJECT

private:
    bool m_cachedIsSlicingActivated;

    // Internal attributes purely related to how the scene is drawn with GL.
    ShaderHelper *m_depthShader;
    ShaderHelper *m_backgroundShader;
    ShaderHelper *m_surfaceFlatShader;
    ShaderHelper *m_surfaceSmoothShader;
    ShaderHelper *m_surfaceTexturedSmoothShader;
    ShaderHelper *m_surfaceTexturedFlatShader;
    ShaderHelper *m_surfaceGridShader;
    ShaderHelper *m_surfaceSliceFlatShader;
    ShaderHelper *m_surfaceSliceSmoothShader;
    ShaderHelper *m_selectionShader;
    float m_heightNormalizer;
    float m_scaleX;
    float m_scaleY;
    float m_scaleZ;
    GLuint m_depthFrameBuffer;
    GLuint m_selectionFrameBuffer;
    GLuint m_selectionDepthBuffer;
    GLuint m_selectionResultTexture;
    GLfloat m_shadowQualityToShader;
    bool m_flatSupported;
    bool m_selectionActive;
    AbstractRenderItem m_dummyRenderItem;
    GLint m_shadowQualityMultiplier;
    QPoint m_selectedPoint;
    QSurface3DSeries *m_selectedSeries;
    QPoint m_clickedPosition;
    bool m_selectionTexturesDirty;
    GLuint m_noShadowTexture;
    bool m_flipHorizontalGrid;

public:
    explicit Surface3DRenderer(Surface3DController *controller);
    ~Surface3DRenderer();

    void updateData();
    void updateSeries(const QList<QAbstract3DSeries *> &seriesList);
    void updateSurfaceTextures(QVector<QSurface3DSeries *> seriesList);
    SeriesRenderCache *createNewCache(QAbstract3DSeries *series);
    void cleanCache(SeriesRenderCache *cache);
    void updateSelectionMode(QAbstract3DGraph::SelectionFlags mode);
    void updateRows(const QVector<Surface3DController::ChangeRow> &rows);
    void updateItems(const QVector<Surface3DController::ChangeItem> &points);
    void updateScene(Q3DScene *scene);
    void updateSlicingActive(bool isSlicing);
    void updateSelectedPoint(const QPoint &position, QSurface3DSeries *series);
    void updateFlipHorizontalGrid(bool flip);
    inline QPoint clickedPosition() const { return m_clickedPosition; }
    void resetClickedStatus();
    QVector3D convertPositionToTranslation(const QVector3D &position, bool isAbsolute);
    void updateAxisLabels(QAbstract3DAxis::AxisOrientation orientation,
                          const QStringList &labels);
    void updateAxisTitleVisibility(QAbstract3DAxis::AxisOrientation orientation,
                                   bool visible);
    void updateMargin(float margin);

    void render(GLuint defaultFboHandle = 0);

protected:
    void contextCleanup();
    void initializeOpenGL();
    virtual void fixCameraTarget(QVector3D &target);
    virtual void getVisibleItemBounds(QVector3D &minBounds, QVector3D &maxBounds);

Q_SIGNALS:
    void flatShadingSupportedChanged(bool supported);

private:
    void checkFlatSupport(SurfaceSeriesRenderCache *cache);
    void updateObjects(SurfaceSeriesRenderCache *cache, bool dimensionChanged);
    void updateSliceDataModel(const QPoint &point);
    QPoint mapCoordsToSampleSpace(SurfaceSeriesRenderCache *cache, const QPointF &coords);
    void findMatchingRow(float z, int &sample, int direction, QSurfaceDataArray &dataArray);
    void findMatchingColumn(float x, int &sample, int direction, QSurfaceDataArray &dataArray);
    void updateSliceObject(SurfaceSeriesRenderCache *cache, const QPoint &point);
    void updateShadowQuality(QAbstract3DGraph::ShadowQuality quality);
    void updateTextures();
    void initShaders(const QString &vertexShader, const QString &fragmentShader);
    QRect calculateSampleRect(const QSurfaceDataArray &array);
    void loadBackgroundMesh();

    void drawSlicedScene();
    void drawScene(GLuint defaultFboHandle);
    void drawLabels(bool drawSelection, const Q3DCamera *activeCamera,
                    const QMatrix4x4 &viewMatrix, const QMatrix4x4 &projectionMatrix);

    void calculateSceneScalingFactors();
    void initBackgroundShaders(const QString &vertexShader, const QString &fragmentShader);
    void initSelectionShaders();
    void initSurfaceShaders();
    void initSelectionBuffer();
    void initDepthShader();
    void updateSelectionTextures();
    void createSelectionTexture(SurfaceSeriesRenderCache *cache, uint &lastSelectionId);
    void idToRGBA(uint id, uchar *r, uchar *g, uchar *b, uchar *a);
    void fillIdCorner(uchar *p, uchar r, uchar g, uchar b, uchar a);
    void surfacePointSelected(const QPoint &point);
    void updateSelectionPoint(SurfaceSeriesRenderCache *cache, const QPoint &point, bool label);
    QPoint selectionIdToSurfacePoint(uint id);
    void updateDepthBuffer();
    void emitSelectedPointChanged(QPoint position);

    Q_DISABLE_COPY(Surface3DRenderer)

    friend class SurfaceObject;
};

QT_END_NAMESPACE_DATAVISUALIZATION

#endif
