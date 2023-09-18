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

#ifndef Q3DBARSRENDERER_P_H
#define Q3DBARSRENDERER_P_H

#include "datavisualizationglobal_p.h"
#include "bars3dcontroller_p.h"
#include "abstract3drenderer_p.h"
#include "barrenderitem_p.h"

QT_BEGIN_NAMESPACE
class QPoint;
class QSizeF;
QT_END_NAMESPACE

QT_BEGIN_NAMESPACE_DATAVISUALIZATION

class ShaderHelper;
class LabelItem;
class Q3DScene;
class BarSeriesRenderCache;

class QT_DATAVISUALIZATION_EXPORT Bars3DRenderer : public Abstract3DRenderer
{
    Q_OBJECT

private:
    // Cached state based on emitted signals from the controller
    QSizeF m_cachedBarThickness;
    QSizeF m_cachedBarSpacing;
    bool m_cachedIsSlicingActivated;
    int m_cachedRowCount;
    int m_cachedColumnCount;

    // Internal state
    BarRenderItem *m_selectedBar; // points to renderitem array
    AxisRenderCache *m_sliceCache; // not owned
    const LabelItem *m_sliceTitleItem; // not owned
    bool m_updateLabels;
    ShaderHelper *m_barShader;
    ShaderHelper *m_barGradientShader;
    ShaderHelper *m_depthShader;
    ShaderHelper *m_selectionShader;
    ShaderHelper *m_backgroundShader;
    GLuint m_bgrTexture;
    GLuint m_selectionTexture;
    GLuint m_depthFrameBuffer;
    GLuint m_selectionFrameBuffer;
    GLuint m_selectionDepthBuffer;
    GLfloat m_shadowQualityToShader;
    GLint m_shadowQualityMultiplier;
    GLfloat m_heightNormalizer;
    GLfloat m_gradientFraction;
    GLfloat m_backgroundAdjustment;
    GLfloat m_rowWidth;
    GLfloat m_columnDepth;
    GLfloat m_maxDimension;
    GLfloat m_scaleX;
    GLfloat m_scaleZ;
    GLfloat m_scaleFactor;
    GLfloat m_maxSceneSize;
    QPoint m_visualSelectedBarPos;
    QPoint m_selectedBarPos;
    BarSeriesRenderCache *m_selectedSeriesCache;
    BarRenderItem m_dummyBarRenderItem;
    bool m_noZeroInRange;
    float m_seriesScaleX;
    float m_seriesScaleZ;
    float m_seriesStep;
    float m_seriesStart;
    QPoint m_clickedPosition;
    bool m_keepSeriesUniform;
    bool m_haveUniformColorSeries;
    bool m_haveGradientSeries;
    float m_zeroPosition;
    float m_xScaleFactor;
    float m_zScaleFactor;
    float m_floorLevel;
    float m_actualFloorLevel;

public:
    explicit Bars3DRenderer(Bars3DController *controller);
    ~Bars3DRenderer();

    void updateData();
    void updateSeries(const QList<QAbstract3DSeries *> &seriesList);
    SeriesRenderCache *createNewCache(QAbstract3DSeries *series);
    void updateRows(const QVector<Bars3DController::ChangeRow> &rows);
    void updateItems(const QVector<Bars3DController::ChangeItem> &items);
    void updateScene(Q3DScene *scene);
    void render(GLuint defaultFboHandle = 0);

    QVector3D convertPositionToTranslation(const QVector3D &position, bool isAbsolute);

    void updateAspectRatio(float ratio);
    void updateFloorLevel(float level);
    void updateMargin(float margin);

protected:
    void contextCleanup();
    virtual void initializeOpenGL();
    virtual void fixCameraTarget(QVector3D &target);
    virtual void getVisibleItemBounds(QVector3D &minBounds, QVector3D &maxBounds);

public Q_SLOTS:
    void updateMultiSeriesScaling(bool uniform);
    void updateBarSpecs(GLfloat thicknessRatio = 1.0f,
                        const QSizeF &spacing = QSizeF(1.0, 1.0),
                        bool relative = true);
    void updateSlicingActive(bool isSlicing);
    void updateSelectedBar(const QPoint &position, QBar3DSeries *series);
    inline QPoint clickedPosition() const { return m_clickedPosition; }
    void resetClickedStatus();

    // Overloaded from abstract renderer
    virtual void updateAxisRange(QAbstract3DAxis::AxisOrientation orientation, float min,
                                 float max);
    virtual void updateAxisReversed(QAbstract3DAxis::AxisOrientation orientation,
                                    bool enable);

private:
    virtual void initShaders(const QString &vertexShader, const QString &fragmentShader);
    virtual void initGradientShaders(const QString &vertexShader, const QString &fragmentShader);
    virtual void updateShadowQuality(QAbstract3DGraph::ShadowQuality quality);
    virtual void updateTextures();
    virtual void fixMeshFileName(QString &fileName, QAbstract3DSeries::Mesh mesh);

    void drawSlicedScene();
    void drawScene(GLuint defaultFboHandle);
    void drawLabels(bool drawSelection, const Q3DCamera *activeCamera,
                    const QMatrix4x4 &viewMatrix, const QMatrix4x4 &projectionMatrix);

    bool drawBars(BarRenderItem **selectedBar, const QMatrix4x4 &depthProjectionViewMatrix,
                  const QMatrix4x4 &projectionViewMatrix, const QMatrix4x4 &viewMatrix,
                  GLint startRow, GLint stopRow, GLint stepRow,
                  GLint startBar, GLint stopBar, GLint stepBar, GLfloat reflection = 1.0f);
    void drawBackground(GLfloat backgroundRotation, const QMatrix4x4 &depthProjectionViewMatrix,
                        const QMatrix4x4 &projectionViewMatrix, const QMatrix4x4 &viewMatrix,
                        bool reflectingDraw = false, bool drawingSelectionBuffer = false);
    void drawGridLines(const QMatrix4x4 &depthProjectionViewMatrix,
                       const QMatrix4x4 &projectionViewMatrix,
                       const QMatrix4x4 &viewMatrix);

    void loadBackgroundMesh();
    void initSelectionShader();
    void initBackgroundShaders(const QString &vertexShader, const QString &fragmentShader);
    void initSelectionBuffer();
    void initDepthShader();
    void updateDepthBuffer();
    void calculateSceneScalingFactors();
    void calculateHeightAdjustment();
    Abstract3DController::SelectionType isSelected(int row, int bar,
                                                   const BarSeriesRenderCache *cache);
    QPoint selectionColorToArrayPosition(const QVector4D &selectionColor);
    QBar3DSeries *selectionColorToSeries(const QVector4D &selectionColor);

    inline void updateRenderRow(const QBarDataRow *dataRow, BarRenderItemRow &renderRow);
    inline void updateRenderItem(const QBarDataItem &dataItem, BarRenderItem &renderItem);

    Q_DISABLE_COPY(Bars3DRenderer)
};

QT_END_NAMESPACE_DATAVISUALIZATION

#endif
