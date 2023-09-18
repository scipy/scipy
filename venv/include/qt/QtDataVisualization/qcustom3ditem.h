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

#ifndef QCUSTOM3DITEM_H
#define QCUSTOM3DITEM_H

#include <QtDataVisualization/qdatavisualizationglobal.h>
#include <QtGui/QImage>
#include <QtGui/QVector3D>
#include <QtGui/QQuaternion>

QT_BEGIN_NAMESPACE_DATAVISUALIZATION

class QCustom3DItemPrivate;

class QT_DATAVISUALIZATION_EXPORT QCustom3DItem : public QObject
{
    Q_OBJECT
    Q_PROPERTY(QString meshFile READ meshFile WRITE setMeshFile NOTIFY meshFileChanged)
    Q_PROPERTY(QString textureFile READ textureFile WRITE setTextureFile NOTIFY textureFileChanged)
    Q_PROPERTY(QVector3D position READ position WRITE setPosition NOTIFY positionChanged)
    Q_PROPERTY(bool positionAbsolute READ isPositionAbsolute WRITE setPositionAbsolute NOTIFY positionAbsoluteChanged)
    Q_PROPERTY(QVector3D scaling READ scaling WRITE setScaling NOTIFY scalingChanged)
    Q_PROPERTY(QQuaternion rotation READ rotation WRITE setRotation NOTIFY rotationChanged)
    Q_PROPERTY(bool visible READ isVisible WRITE setVisible NOTIFY visibleChanged)
    Q_PROPERTY(bool shadowCasting READ isShadowCasting WRITE setShadowCasting NOTIFY shadowCastingChanged)
    Q_PROPERTY(bool scalingAbsolute READ isScalingAbsolute WRITE setScalingAbsolute NOTIFY scalingAbsoluteChanged REVISION 1)

public:
    explicit QCustom3DItem(QObject *parent = nullptr);
    explicit QCustom3DItem(const QString &meshFile, const QVector3D &position,
                           const QVector3D &scaling, const QQuaternion &rotation,
                           const QImage &texture, QObject *parent = nullptr);
    virtual ~QCustom3DItem();

    void setMeshFile(const QString &meshFile);
    QString meshFile() const;

    void setTextureFile(const QString &textureFile);
    QString textureFile() const;

    void setPosition(const QVector3D &position);
    QVector3D position() const;

    void setPositionAbsolute(bool positionAbsolute);
    bool isPositionAbsolute() const;

    void setScaling(const QVector3D &scaling);
    QVector3D scaling() const;

    void setScalingAbsolute(bool scalingAbsolute);
    bool isScalingAbsolute() const;

    void setRotation(const QQuaternion &rotation);
    QQuaternion rotation();

    void setVisible(bool visible);
    bool isVisible() const;

    void setShadowCasting(bool enabled);
    bool isShadowCasting() const;

    Q_INVOKABLE void setRotationAxisAndAngle(const QVector3D &axis, float angle);

    void setTextureImage(const QImage &textureImage);

Q_SIGNALS:
    void meshFileChanged(const QString &meshFile);
    void textureFileChanged(const QString &textureFile);
    void positionChanged(const QVector3D &position);
    void positionAbsoluteChanged(bool positionAbsolute);
    void scalingChanged(const QVector3D &scaling);
    void rotationChanged(const QQuaternion &rotation);
    void visibleChanged(bool visible);
    void shadowCastingChanged(bool shadowCasting);
    Q_REVISION(1) void scalingAbsoluteChanged(bool scalingAbsolute);

protected:
    QCustom3DItem(QCustom3DItemPrivate *d, QObject *parent = nullptr);

    QScopedPointer<QCustom3DItemPrivate> d_ptr;

private:
    Q_DISABLE_COPY(QCustom3DItem)

    friend class Abstract3DRenderer;
    friend class Abstract3DController;
};

QT_END_NAMESPACE_DATAVISUALIZATION

#endif
