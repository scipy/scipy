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

#ifndef ABSTRACT3DRENDERER_P_H
#define ABSTRACT3DRENDERER_P_H

#include <QtGui/QOpenGLFunctions>
#if !defined(QT_OPENGL_ES_2)
#  include <QtGui/QOpenGLFunctions_2_1>
#endif
#include "datavisualizationglobal_p.h"
#include "abstract3dcontroller_p.h"
#include "axisrendercache_p.h"
#include "seriesrendercache_p.h"
#include "customrenderitem_p.h"

QT_FORWARD_DECLARE_CLASS(QOffscreenSurface)

QT_BEGIN_NAMESPACE_DATAVISUALIZATION

class TextureHelper;
class Theme;
class Drawer;

class Abstract3DRenderer : public QObject, protected QOpenGLFunctions
{
    Q_OBJECT

protected:
    enum SelectionState {
        SelectNone = 0,
        SelectOnScene,
        SelectOnOverview,
        SelectOnSlice
    };

    enum RenderingState {
        RenderingNormal = 0,
        RenderingSelection,
        RenderingDepth
    };

public:
    virtual ~Abstract3DRenderer();

    virtual void updateData() = 0;
    virtual void updateSeries(const QList<QAbstract3DSeries *> &seriesList);
    virtual void updateCustomData(const QList<QCustom3DItem *> &customItems);
    virtual void updateCustomItems();
    virtual void updateCustomItemPositions();
    virtual SeriesRenderCache *createNewCache(QAbstract3DSeries *series);
    virtual void cleanCache(SeriesRenderCache *cache);
    virtual void render(GLuint defaultFboHandle);

    virtual void updateTheme(Q3DTheme *theme);
    virtual void updateSelectionMode(QAbstract3DGraph::SelectionFlags newMode);
    virtual void updateOptimizationHint(QAbstract3DGraph::OptimizationHints hint);
    virtual void updateScene(Q3DScene *scene);
    virtual void updateTextures();
    virtual void initSelectionBuffer() = 0;
    virtual void updateSelectionState(SelectionState state);

    virtual void updateDepthBuffer() = 0;
    virtual void updateShadowQuality(QAbstract3DGraph::ShadowQuality quality) = 0;
    virtual void initShaders(const QString &vertexShader, const QString &fragmentShader) = 0;
    virtual void initGradientShaders(const QString &vertexShader, const QString &fragmentShader);
    virtual void initStaticSelectedItemShaders(const QString &vertexShader,
                                               const QString &fragmentShader,
                                               const QString &gradientVertexShader,
                                               const QString &gradientFragmentShader);
    virtual void initBackgroundShaders(const QString &vertexShader,
                                       const QString &fragmentShader) = 0;
    virtual void initCustomItemShaders(const QString &vertexShader,
                                       const QString &fragmentShader);
    virtual void initVolumeTextureShaders(const QString &vertexShader,
                                          const QString &fragmentShader,
                                          const QString &fragmentLowDefShader,
                                          const QString &sliceShader,
                                          const QString &sliceFrameVertexShader,
                                          const QString &sliceFrameShader);
    virtual void initLabelShaders(const QString &vertexShader, const QString &fragmentShader);
    virtual void initCursorPositionShaders(const QString &vertexShader,
                                           const QString &fragmentShader);
    virtual void initCursorPositionBuffer();

    virtual void updateAxisType(QAbstract3DAxis::AxisOrientation orientation,
                                QAbstract3DAxis::AxisType type);
    virtual void updateAxisTitle(QAbstract3DAxis::AxisOrientation orientation,
                                 const QString &title);
    virtual void updateAxisLabels(QAbstract3DAxis::AxisOrientation orientation,
                                  const QStringList &labels);
    virtual void updateAxisRange(QAbstract3DAxis::AxisOrientation orientation, float min,
                                 float max);
    virtual void updateAxisSegmentCount(QAbstract3DAxis::AxisOrientation orientation, int count);
    virtual void updateAxisSubSegmentCount(QAbstract3DAxis::AxisOrientation orientation,
                                           int count);
    virtual void updateAxisLabelFormat(QAbstract3DAxis::AxisOrientation orientation,
                                       const QString &format);
    virtual void updateAxisReversed(QAbstract3DAxis::AxisOrientation orientation,
                                    bool enable);
    virtual void updateAxisFormatter(QAbstract3DAxis::AxisOrientation orientation,
                                     QValue3DAxisFormatter *formatter);
    virtual void updateAxisLabelAutoRotation(QAbstract3DAxis::AxisOrientation orientation,
                                             float angle);
    virtual void updateAxisTitleVisibility(QAbstract3DAxis::AxisOrientation orientation,
                                           bool visible);
    virtual void updateAxisTitleFixed(QAbstract3DAxis::AxisOrientation orientation,
                                      bool fixed);
    virtual void modifiedSeriesList(const QVector<QAbstract3DSeries *> &seriesList);

    virtual void fixMeshFileName(QString &fileName, QAbstract3DSeries::Mesh mesh);

    virtual CustomRenderItem *addCustomItem(QCustom3DItem *item);
    virtual void updateCustomItem(CustomRenderItem *renderItem);

    virtual void updateAspectRatio(float ratio);
    virtual void updateHorizontalAspectRatio(float ratio);
    virtual void updatePolar(bool enable);
    virtual void updateRadialLabelOffset(float offset);
    virtual void updateMargin(float margin);

    virtual QVector3D convertPositionToTranslation(const QVector3D &position,
                                                   bool isAbsolute) = 0;

    void generateBaseColorTexture(const QColor &color, GLuint *texture);
    void fixGradientAndGenerateTexture(QLinearGradient *gradient, GLuint *gradientTexture);

    inline bool isClickQueryResolved() const { return m_clickResolved; }
    inline void clearClickQueryResolved() { m_clickResolved = false; }
    inline QPoint cachedClickQuery() const { return m_cachedScene->selectionQueryPosition(); }
    inline QAbstract3DSeries *clickedSeries() const { return m_clickedSeries; }
    inline QAbstract3DGraph::ElementType clickedType() { return m_clickedType; }
    inline bool isGraphPositionQueryResolved() const { return m_graphPositionQueryResolved; }
    inline void clearGraphPositionQueryResolved() { m_graphPositionQueryResolved = false; }
    inline QVector3D queriedGraphPosition() const { return m_queriedGraphPosition; }
    inline QPoint cachedGraphPositionQuery() const { return m_cachedScene->graphPositionQuery(); }

    LabelItem &selectionLabelItem();
    void setSelectionLabel(const QString &label);
    QString &selectionLabel();

    void drawCustomItems(RenderingState state, ShaderHelper *regularShader,
                         const QMatrix4x4 &viewMatrix,
                         const QMatrix4x4 &projectionViewMatrix,
                         const QMatrix4x4 &depthProjectionViewMatrix,
                         GLuint depthTexture, GLfloat shadowQuality, GLfloat reflection = 1.0f);

    QVector4D indexToSelectionColor(GLint index);
    void calculatePolarXZ(const QVector3D &dataPos, float &x, float &z) const;

Q_SIGNALS:
    void needRender(); // Emit this if something in renderer causes need for another render pass.
    void requestShadowQuality(QAbstract3DGraph::ShadowQuality quality); // For automatic quality adjustments

protected:
    Abstract3DRenderer(Abstract3DController *controller);

    virtual void contextCleanup();
    virtual void initializeOpenGL();

    void reInitShaders();
    virtual void handleShadowQualityChange();
    virtual void handleResize();

    AxisRenderCache &axisCacheForOrientation(QAbstract3DAxis::AxisOrientation orientation);

    virtual void lowerShadowQuality();

    void fixGradient(QLinearGradient *gradient, GLuint *gradientTexture);

    void calculateZoomLevel();
    void drawAxisTitleY(const QVector3D &sideLabelRotation, const QVector3D &backLabelRotation,
                        const QVector3D &sideLabelTrans, const QVector3D &backLabelTrans,
                        const QQuaternion &totalSideRotation, const QQuaternion &totalBackRotation,
                        AbstractRenderItem &dummyItem, const Q3DCamera *activeCamera,
                        float labelsMaxWidth,
                        const QMatrix4x4 &viewMatrix, const QMatrix4x4 &projectionMatrix,
                        ShaderHelper *shader);
    void drawAxisTitleX(const QVector3D &labelRotation, const QVector3D &labelTrans,
                        const QQuaternion &totalRotation, AbstractRenderItem &dummyItem,
                        const Q3DCamera *activeCamera, float labelsMaxWidth,
                        const QMatrix4x4 &viewMatrix, const QMatrix4x4 &projectionMatrix,
                        ShaderHelper *shader, bool radial = false);
    void drawAxisTitleZ(const QVector3D &labelRotation, const QVector3D &labelTrans,
                        const QQuaternion &totalRotation, AbstractRenderItem &dummyItem,
                        const Q3DCamera *activeCamera, float labelsMaxWidth,
                        const QMatrix4x4 &viewMatrix, const QMatrix4x4 &projectionMatrix,
                        ShaderHelper *shader);

    void loadGridLineMesh();
    void loadLabelMesh();
    void loadPositionMapperMesh();

    void drawRadialGrid(ShaderHelper *shader, float yFloorLinePos,
                        const QMatrix4x4 &projectionViewMatrix, const QMatrix4x4 &depthMatrix);
    void drawAngularGrid(ShaderHelper *shader, float yFloorLinePos,
                         const QMatrix4x4 &projectionViewMatrix, const QMatrix4x4 &depthMatrix);

    float calculatePolarBackgroundMargin();
    virtual void fixCameraTarget(QVector3D &target) = 0;
    void updateCameraViewport();

    void recalculateCustomItemScalingAndPos(CustomRenderItem *item);
    virtual void getVisibleItemBounds(QVector3D &minBounds, QVector3D &maxBounds) = 0;
    void drawVolumeSliceFrame(const CustomRenderItem *item, Qt::Axis axis,
                              const QMatrix4x4 &projectionViewMatrix);
    void queriedGraphPosition(const QMatrix4x4 &projectionViewMatrix, const QVector3D &scaling,
                              GLuint defaultFboHandle);

    bool m_hasNegativeValues;
    Q3DTheme *m_cachedTheme;
    Drawer *m_drawer;
    QRect m_viewport;
    QAbstract3DGraph::ShadowQuality m_cachedShadowQuality;
    GLfloat m_autoScaleAdjustment;

    QAbstract3DGraph::SelectionFlags m_cachedSelectionMode;
    QAbstract3DGraph::OptimizationHints m_cachedOptimizationHint;

    AxisRenderCache m_axisCacheX;
    AxisRenderCache m_axisCacheY;
    AxisRenderCache m_axisCacheZ;
    TextureHelper *m_textureHelper;
    GLuint m_depthTexture;

    Q3DScene *m_cachedScene;
    bool m_selectionDirty;
    SelectionState m_selectionState;
    QPoint m_inputPosition;
    QHash<QAbstract3DSeries *, SeriesRenderCache *> m_renderCacheList;
    CustomRenderItemArray m_customRenderCache;
    QList<QCustom3DItem *> m_customItemDrawOrder;
    QRect m_primarySubViewport;
    QRect m_secondarySubViewport;
    float m_devicePixelRatio;
    bool m_selectionLabelDirty;
    bool m_clickResolved;
    bool m_graphPositionQueryPending;
    bool m_graphPositionQueryResolved;
    QAbstract3DSeries *m_clickedSeries;
    QAbstract3DGraph::ElementType m_clickedType;
    int m_selectedLabelIndex;
    int m_selectedCustomItemIndex;
    QVector3D m_queriedGraphPosition;
    QPoint m_graphPositionQuery;

    QString m_selectionLabel;
    LabelItem *m_selectionLabelItem;
    int m_visibleSeriesCount;

    ShaderHelper *m_customItemShader;
    ShaderHelper *m_volumeTextureShader;
    ShaderHelper *m_volumeTextureLowDefShader;
    ShaderHelper *m_volumeTextureSliceShader;
    ShaderHelper *m_volumeSliceFrameShader;
    ShaderHelper *m_labelShader;
    ShaderHelper *m_cursorPositionShader;
    GLuint m_cursorPositionFrameBuffer;
    GLuint m_cursorPositionTexture;

    bool m_useOrthoProjection;
    bool m_xFlipped;
    bool m_yFlipped;
    bool m_zFlipped;
    bool m_yFlippedForGrid;

    ObjectHelper *m_backgroundObj; // Shared reference
    ObjectHelper *m_gridLineObj; // Shared reference
    ObjectHelper *m_labelObj; // Shared reference
    ObjectHelper *m_positionMapperObj; // Shared reference

    float m_graphAspectRatio;
    float m_graphHorizontalAspectRatio;
    bool m_polarGraph;
    float m_radialLabelOffset;
    float m_polarRadius;

    QQuaternion m_xRightAngleRotation;
    QQuaternion m_yRightAngleRotation;
    QQuaternion m_zRightAngleRotation;
    QQuaternion m_xRightAngleRotationNeg;
    QQuaternion m_yRightAngleRotationNeg;
    QQuaternion m_zRightAngleRotationNeg;
    QQuaternion m_xFlipRotation;
    QQuaternion m_zFlipRotation;

    float m_requestedMargin;
    float m_vBackgroundMargin;
    float m_hBackgroundMargin;
    float m_scaleXWithBackground;
    float m_scaleYWithBackground;
    float m_scaleZWithBackground;

    QVector3D m_oldCameraTarget;

    bool m_reflectionEnabled;
    qreal m_reflectivity;

    QLocale m_locale;
#if !defined(QT_OPENGL_ES_2)
    QOpenGLFunctions_2_1 *m_funcs_2_1;  // Not owned
#endif
    QPointer<QOpenGLContext> m_context; // Not owned
    bool m_isOpenGLES;

private:
    friend class Abstract3DController;
};

QT_END_NAMESPACE_DATAVISUALIZATION

#endif
