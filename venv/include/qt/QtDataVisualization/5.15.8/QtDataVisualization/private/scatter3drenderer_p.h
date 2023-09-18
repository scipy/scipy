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

#ifndef Q3DSCATTERRENDERER_P_H
#define Q3DSCATTERRENDERER_P_H

#include "datavisualizationglobal_p.h"
#include "scatter3dcontroller_p.h"
#include "abstract3drenderer_p.h"
#include "scatterrenderitem_p.h"

QT_FORWARD_DECLARE_CLASS(QSizeF)

QT_BEGIN_NAMESPACE_DATAVISUALIZATION

class ShaderHelper;
class Q3DScene;
class ScatterSeriesRenderCache;
class QScatterDataItem;

class QT_DATAVISUALIZATION_EXPORT Scatter3DRenderer : public Abstract3DRenderer
{
    Q_OBJECT

private:
    // Internal state
    ScatterRenderItem *m_selectedItem; // points to renderitem array
    bool m_updateLabels;
    ShaderHelper *m_dotShader;
    ShaderHelper *m_dotGradientShader;
    ShaderHelper *m_staticSelectedItemGradientShader;
    ShaderHelper *m_staticSelectedItemShader;
    ShaderHelper *m_pointShader;
    ShaderHelper *m_depthShader;
    ShaderHelper *m_selectionShader;
    ShaderHelper *m_backgroundShader;
    ShaderHelper *m_staticGradientPointShader;
    GLuint m_bgrTexture;
    GLuint m_selectionTexture;
    GLuint m_depthFrameBuffer;
    GLuint m_selectionFrameBuffer;
    GLuint m_selectionDepthBuffer;
    GLfloat m_shadowQualityToShader;
    GLint m_shadowQualityMultiplier;
    float m_scaleX;
    float m_scaleY;
    float m_scaleZ;
    int m_selectedItemIndex;
    ScatterSeriesRenderCache *m_selectedSeriesCache;
    ScatterSeriesRenderCache *m_oldSelectedSeriesCache;
    GLfloat m_dotSizeScale;
    ScatterRenderItem m_dummyRenderItem;
    GLfloat m_maxItemSize;
    int m_clickedIndex;
    bool m_havePointSeries;
    bool m_haveMeshSeries;
    bool m_haveUniformColorMeshSeries;
    bool m_haveGradientMeshSeries;

public:
    explicit Scatter3DRenderer(Scatter3DController *controller);
    ~Scatter3DRenderer();

    void updateData();
    void updateSeries(const QList<QAbstract3DSeries *> &seriesList);
    SeriesRenderCache *createNewCache(QAbstract3DSeries *series);
    void updateItems(const QVector<Scatter3DController::ChangeItem> &items);
    void updateScene(Q3DScene *scene);
    void updateAxisLabels(QAbstract3DAxis::AxisOrientation orientation,
                          const QStringList &labels);
    void updateAxisTitleVisibility(QAbstract3DAxis::AxisOrientation orientation,
                                   bool visible);
    void updateOptimizationHint(QAbstract3DGraph::OptimizationHints hint);
    void updateMargin(float margin);

    QVector3D convertPositionToTranslation(const QVector3D &position, bool isAbsolute);

    inline int clickedIndex() const { return m_clickedIndex; }
    void resetClickedStatus();

    void render(GLuint defaultFboHandle);

public Q_SLOTS:
    void updateSelectedItem(int index, QScatter3DSeries *series);

protected:
    void contextCleanup();
    virtual void initializeOpenGL();
    virtual void fixCameraTarget(QVector3D &target);
    virtual void getVisibleItemBounds(QVector3D &minBounds, QVector3D &maxBounds);

private:
    virtual void initShaders(const QString &vertexShader, const QString &fragmentShader);
    virtual void initGradientShaders(const QString &vertexShader, const QString &fragmentShader);
    virtual void initStaticSelectedItemShaders(const QString &vertexShader,
                                               const QString &fragmentShader,
                                               const QString &gradientVertexShader,
                                               const QString &gradientFragmentShader);
    virtual void updateShadowQuality(QAbstract3DGraph::ShadowQuality quality);
    virtual void updateTextures();
    virtual void fixMeshFileName(QString &fileName, QAbstract3DSeries::Mesh mesh);

    void drawScene(GLuint defaultFboHandle);
    void drawLabels(bool drawSelection, const Q3DCamera *activeCamera,
                    const QMatrix4x4 &viewMatrix, const QMatrix4x4 &projectionMatrix);

    void loadBackgroundMesh();
    void initSelectionShader();
    void initBackgroundShaders(const QString &vertexShader, const QString &fragmentShader);
    void initStaticPointShaders(const QString &vertexShader, const QString &fragmentShader);
    void initSelectionBuffer();
    void initDepthShader();
    void updateDepthBuffer();
    void initPointShader();
    void calculateTranslation(ScatterRenderItem &item);
    void calculateSceneScalingFactors();

    void selectionColorToSeriesAndIndex(const QVector4D &color, int &index,
                                        QAbstract3DSeries *&series);
    inline void updateRenderItem(const QScatterDataItem &dataItem, ScatterRenderItem &renderItem);

    Q_DISABLE_COPY(Scatter3DRenderer)
};

QT_END_NAMESPACE_DATAVISUALIZATION

#endif
