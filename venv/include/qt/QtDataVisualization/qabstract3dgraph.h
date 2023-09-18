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

#ifndef QABSTRACT3DGRAPH_H
#define QABSTRACT3DGRAPH_H

#include <QtDataVisualization/qdatavisualizationglobal.h>
#include <QtDataVisualization/q3dtheme.h>
#include <QtDataVisualization/q3dscene.h>
#include <QtDataVisualization/qabstract3dinputhandler.h>
#include <QtGui/QWindow>
#include <QtGui/QOpenGLFunctions>
#include <QtCore/QLocale>

QT_BEGIN_NAMESPACE_DATAVISUALIZATION

class QAbstract3DGraphPrivate;
class QCustom3DItem;
class QAbstract3DAxis;

class QT_DATAVISUALIZATION_EXPORT QAbstract3DGraph : public QWindow, protected QOpenGLFunctions
{
    Q_OBJECT
    Q_ENUMS(ShadowQuality)
    Q_ENUMS(ElementType)
    Q_FLAGS(SelectionFlag SelectionFlags)
    Q_FLAGS(OptimizationHint OptimizationHints)
    Q_PROPERTY(QAbstract3DInputHandler* activeInputHandler READ activeInputHandler WRITE setActiveInputHandler NOTIFY activeInputHandlerChanged)
    Q_PROPERTY(Q3DTheme* activeTheme READ activeTheme WRITE setActiveTheme NOTIFY activeThemeChanged)
    Q_PROPERTY(SelectionFlags selectionMode READ selectionMode WRITE setSelectionMode NOTIFY selectionModeChanged)
    Q_PROPERTY(ShadowQuality shadowQuality READ shadowQuality WRITE setShadowQuality NOTIFY shadowQualityChanged)
    Q_PROPERTY(Q3DScene* scene READ scene)
    Q_PROPERTY(bool measureFps READ measureFps WRITE setMeasureFps NOTIFY measureFpsChanged)
    Q_PROPERTY(qreal currentFps READ currentFps NOTIFY currentFpsChanged)
    Q_PROPERTY(bool orthoProjection READ isOrthoProjection WRITE setOrthoProjection NOTIFY orthoProjectionChanged)
    Q_PROPERTY(ElementType selectedElement READ selectedElement NOTIFY selectedElementChanged)
    Q_PROPERTY(qreal aspectRatio READ aspectRatio WRITE setAspectRatio NOTIFY aspectRatioChanged)
    Q_PROPERTY(OptimizationHints optimizationHints READ optimizationHints WRITE setOptimizationHints NOTIFY optimizationHintsChanged)
    Q_PROPERTY(bool polar READ isPolar WRITE setPolar NOTIFY polarChanged)
    Q_PROPERTY(float radialLabelOffset READ radialLabelOffset WRITE setRadialLabelOffset NOTIFY radialLabelOffsetChanged)
    Q_PROPERTY(qreal horizontalAspectRatio READ horizontalAspectRatio WRITE setHorizontalAspectRatio NOTIFY horizontalAspectRatioChanged)
    Q_PROPERTY(bool reflection READ isReflection WRITE setReflection NOTIFY reflectionChanged)
    Q_PROPERTY(qreal reflectivity READ reflectivity WRITE setReflectivity NOTIFY reflectivityChanged)
    Q_PROPERTY(QLocale locale READ locale WRITE setLocale NOTIFY localeChanged)
    Q_PROPERTY(QVector3D queriedGraphPosition READ queriedGraphPosition NOTIFY queriedGraphPositionChanged)
    Q_PROPERTY(qreal margin READ margin WRITE setMargin NOTIFY marginChanged)

protected:
    explicit QAbstract3DGraph(QAbstract3DGraphPrivate *d, const QSurfaceFormat *format,
                              QWindow *parent = nullptr);

public:
    enum SelectionFlag {
        SelectionNone              = 0,
        SelectionItem              = 1,
        SelectionRow               = 2,
        SelectionItemAndRow        = SelectionItem | SelectionRow,
        SelectionColumn            = 4,
        SelectionItemAndColumn     = SelectionItem | SelectionColumn,
        SelectionRowAndColumn      = SelectionRow | SelectionColumn,
        SelectionItemRowAndColumn  = SelectionItem | SelectionRow | SelectionColumn,
        SelectionSlice             = 8,
        SelectionMultiSeries       = 16
    };
    Q_DECLARE_FLAGS(SelectionFlags, SelectionFlag)

    enum ShadowQuality {
        ShadowQualityNone = 0,
        ShadowQualityLow,
        ShadowQualityMedium,
        ShadowQualityHigh,
        ShadowQualitySoftLow,
        ShadowQualitySoftMedium,
        ShadowQualitySoftHigh
    };

    enum ElementType {
        ElementNone = 0,
        ElementSeries,
        ElementAxisXLabel,
        ElementAxisYLabel,
        ElementAxisZLabel,
        ElementCustomItem
    };

    enum OptimizationHint {
        OptimizationDefault = 0,
        OptimizationStatic  = 1
    };
    Q_DECLARE_FLAGS(OptimizationHints, OptimizationHint)

public:
    virtual ~QAbstract3DGraph();

    void addInputHandler(QAbstract3DInputHandler *inputHandler);
    void releaseInputHandler(QAbstract3DInputHandler *inputHandler);
    void setActiveInputHandler(QAbstract3DInputHandler *inputHandler);
    QAbstract3DInputHandler *activeInputHandler() const;
    QList<QAbstract3DInputHandler *> inputHandlers() const;

    void addTheme(Q3DTheme *theme);
    void releaseTheme(Q3DTheme *theme);
    void setActiveTheme(Q3DTheme *theme);
    Q3DTheme *activeTheme() const;
    QList<Q3DTheme *> themes() const;

    void setSelectionMode(SelectionFlags mode);
    SelectionFlags selectionMode() const;

    void setShadowQuality(ShadowQuality quality);
    ShadowQuality shadowQuality() const;
    virtual bool shadowsSupported() const;

    Q3DScene *scene() const;

    void clearSelection();

    int addCustomItem(QCustom3DItem *item);
    void removeCustomItems();
    void removeCustomItem(QCustom3DItem *item);
    void removeCustomItemAt(const QVector3D &position);
    void releaseCustomItem(QCustom3DItem *item);
    QList<QCustom3DItem *> customItems() const;

    int selectedLabelIndex() const;
    QAbstract3DAxis *selectedAxis() const;

    int selectedCustomItemIndex() const;
    QCustom3DItem *selectedCustomItem() const;

    QImage renderToImage(int msaaSamples = 0, const QSize &imageSize = QSize());

    void setMeasureFps(bool enable);
    bool measureFps() const;
    qreal currentFps() const;

    void setOrthoProjection(bool enable);
    bool isOrthoProjection() const;

    ElementType selectedElement() const;

    void setAspectRatio(qreal ratio);
    qreal aspectRatio() const;

    void setOptimizationHints(OptimizationHints hints);
    OptimizationHints optimizationHints() const;

    void setPolar(bool enable);
    bool isPolar() const;

    void setRadialLabelOffset(float offset);
    float radialLabelOffset() const;

    void setHorizontalAspectRatio(qreal ratio);
    qreal horizontalAspectRatio() const;

    void setReflection(bool enable);
    bool isReflection() const;

    void setReflectivity(qreal reflectivity);
    qreal reflectivity() const;

    void setLocale(const QLocale &locale);
    QLocale locale() const;

    QVector3D queriedGraphPosition() const;

    void setMargin(qreal margin);
    qreal margin() const;

    bool hasContext() const;

protected:
    bool event(QEvent *event);
    void resizeEvent(QResizeEvent *event);
    void exposeEvent(QExposeEvent *event);

    void mouseDoubleClickEvent(QMouseEvent *event);
    void touchEvent(QTouchEvent *event);
    void mousePressEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
#if QT_CONFIG(wheelevent)
    void wheelEvent(QWheelEvent *event);
#endif

Q_SIGNALS:
    void activeInputHandlerChanged(QAbstract3DInputHandler *inputHandler);
    void activeThemeChanged(Q3DTheme *theme);
    void selectionModeChanged(QAbstract3DGraph::SelectionFlags mode);
    void shadowQualityChanged(QAbstract3DGraph::ShadowQuality quality);
    void selectedElementChanged(QAbstract3DGraph::ElementType type);
    void measureFpsChanged(bool enabled);
    void currentFpsChanged(qreal fps);
    void orthoProjectionChanged(bool enabled);
    void aspectRatioChanged(qreal ratio);
    void optimizationHintsChanged(QAbstract3DGraph::OptimizationHints hints);
    void polarChanged(bool enabled);
    void radialLabelOffsetChanged(float offset);
    void horizontalAspectRatioChanged(qreal ratio);
    void reflectionChanged(bool enabled);
    void reflectivityChanged(qreal reflectivity);
    void localeChanged(const QLocale &locale);
    void queriedGraphPositionChanged(const QVector3D &data);
    void marginChanged(qreal margin);

private:
    Q_DISABLE_COPY(QAbstract3DGraph)
    QScopedPointer<QAbstract3DGraphPrivate> d_ptr;

    friend class Q3DBars;
    friend class Q3DScatter;
    friend class Q3DSurface;
};
Q_DECLARE_OPERATORS_FOR_FLAGS(QAbstract3DGraph::SelectionFlags)
Q_DECLARE_OPERATORS_FOR_FLAGS(QAbstract3DGraph::OptimizationHints)

QT_END_NAMESPACE_DATAVISUALIZATION

#endif
