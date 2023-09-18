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

#ifndef Q3DSCENE_H
#define Q3DSCENE_H

#include <QtDataVisualization/qdatavisualizationglobal.h>
#include <QtDataVisualization/q3dcamera.h>
#include <QtDataVisualization/q3dlight.h>
#include <QtCore/QObject>
#include <QtCore/QRect>

QT_BEGIN_NAMESPACE_DATAVISUALIZATION

class Q3DScenePrivate;

class QT_DATAVISUALIZATION_EXPORT Q3DScene : public QObject
{
    Q_OBJECT
    Q_PROPERTY(QRect viewport READ viewport NOTIFY viewportChanged)
    Q_PROPERTY(QRect primarySubViewport READ primarySubViewport WRITE setPrimarySubViewport NOTIFY primarySubViewportChanged)
    Q_PROPERTY(QRect secondarySubViewport READ secondarySubViewport WRITE setSecondarySubViewport NOTIFY secondarySubViewportChanged)
    Q_PROPERTY(QPoint selectionQueryPosition READ selectionQueryPosition WRITE setSelectionQueryPosition NOTIFY selectionQueryPositionChanged)
    Q_PROPERTY(bool secondarySubviewOnTop READ isSecondarySubviewOnTop  WRITE setSecondarySubviewOnTop  NOTIFY secondarySubviewOnTopChanged)
    Q_PROPERTY(bool slicingActive READ isSlicingActive WRITE setSlicingActive NOTIFY slicingActiveChanged)
    Q_PROPERTY(Q3DCamera* activeCamera READ activeCamera WRITE setActiveCamera NOTIFY activeCameraChanged)
    Q_PROPERTY(Q3DLight* activeLight READ activeLight WRITE setActiveLight NOTIFY activeLightChanged)
    Q_PROPERTY(float devicePixelRatio READ devicePixelRatio WRITE setDevicePixelRatio NOTIFY devicePixelRatioChanged)
    Q_PROPERTY(QPoint graphPositionQuery READ graphPositionQuery WRITE setGraphPositionQuery NOTIFY graphPositionQueryChanged REVISION 1)

public:
    explicit Q3DScene(QObject *parent = nullptr);
    virtual ~Q3DScene();

    QRect viewport() const;

    QRect primarySubViewport() const;
    void setPrimarySubViewport(const QRect &primarySubViewport);
    bool isPointInPrimarySubView(const QPoint &point);

    QRect secondarySubViewport() const;
    void setSecondarySubViewport(const QRect &secondarySubViewport);
    bool isPointInSecondarySubView(const QPoint &point);

    void setSelectionQueryPosition(const QPoint &point);
    QPoint selectionQueryPosition() const;
    static QPoint invalidSelectionPoint();

    void setGraphPositionQuery(const QPoint &point);
    QPoint graphPositionQuery() const;

    void setSlicingActive(bool isSlicing);
    bool isSlicingActive() const;

    void setSecondarySubviewOnTop(bool isSecondaryOnTop);
    bool isSecondarySubviewOnTop() const;

    Q3DCamera *activeCamera() const;
    void setActiveCamera(Q3DCamera *camera);

    Q3DLight *activeLight() const;
    void setActiveLight(Q3DLight *light);

    float devicePixelRatio() const;
    void setDevicePixelRatio(float pixelRatio);

Q_SIGNALS:
    void viewportChanged(const QRect &viewport);
    void primarySubViewportChanged(const QRect &subViewport);
    void secondarySubViewportChanged(const QRect &subViewport);
    void secondarySubviewOnTopChanged(bool isSecondaryOnTop);
    void slicingActiveChanged(bool isSlicingActive);
    void activeCameraChanged(Q3DCamera *camera);
    void activeLightChanged(Q3DLight *light);
    void devicePixelRatioChanged(float pixelRatio);
    void selectionQueryPositionChanged(const QPoint &position);
    Q_REVISION(1) void graphPositionQueryChanged(const QPoint &position);

private:
    QScopedPointer<Q3DScenePrivate> d_ptr;

    Q_DISABLE_COPY(Q3DScene)

    friend class AbstractDeclarative;
    friend class QAbstract3DGraph;
    friend class QAbstract3DGraphPrivate;
    friend class Abstract3DController;
    friend class Bars3DController;
    friend class Q3DScenePrivate;
    friend class Abstract3DRenderer;
    friend class Bars3DRenderer;
    friend class Surface3DRenderer;
    friend class Scatter3DRenderer;
    friend class Q3DCameraPrivate;
    friend class Q3DObject;
};

QT_END_NAMESPACE_DATAVISUALIZATION

#endif
