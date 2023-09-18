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

#ifndef QABSTRACT3DAXIS_H
#define QABSTRACT3DAXIS_H

#include <QtDataVisualization/qdatavisualizationglobal.h>
#include <QtCore/QObject>
#include <QtCore/QScopedPointer>
#include <QtCore/QStringList>

QT_BEGIN_NAMESPACE_DATAVISUALIZATION

class QAbstract3DAxisPrivate;

class QT_DATAVISUALIZATION_EXPORT QAbstract3DAxis : public QObject
{
    Q_OBJECT
    Q_ENUMS(AxisOrientation)
    Q_ENUMS(AxisType)
    Q_PROPERTY(QString title READ title WRITE setTitle NOTIFY titleChanged)
    Q_PROPERTY(QStringList labels READ labels WRITE setLabels NOTIFY labelsChanged)
    Q_PROPERTY(AxisOrientation orientation READ orientation NOTIFY orientationChanged)
    Q_PROPERTY(AxisType type READ type CONSTANT)
    Q_PROPERTY(float min READ min WRITE setMin NOTIFY minChanged)
    Q_PROPERTY(float max READ max WRITE setMax NOTIFY maxChanged)
    Q_PROPERTY(bool autoAdjustRange READ isAutoAdjustRange WRITE setAutoAdjustRange NOTIFY autoAdjustRangeChanged)
    Q_PROPERTY(float labelAutoRotation READ labelAutoRotation WRITE setLabelAutoRotation NOTIFY labelAutoRotationChanged REVISION 1)
    Q_PROPERTY(bool titleVisible READ isTitleVisible WRITE setTitleVisible NOTIFY titleVisibilityChanged REVISION 1)
    Q_PROPERTY(bool titleFixed READ isTitleFixed WRITE setTitleFixed NOTIFY titleFixedChanged REVISION 1)

public:
    enum AxisOrientation {
        AxisOrientationNone = 0,
        AxisOrientationX = 1,
        AxisOrientationY = 2,
        AxisOrientationZ = 4
    };

    enum AxisType {
        AxisTypeNone = 0,
        AxisTypeCategory = 1,
        AxisTypeValue = 2
    };

protected:
    explicit QAbstract3DAxis(QAbstract3DAxisPrivate *d, QObject *parent = nullptr);

public:
    virtual ~QAbstract3DAxis();

    void setTitle(const QString &title);
    QString title() const;

    void setLabels(const QStringList &labels);
    QStringList labels() const;

    AxisOrientation orientation() const;
    AxisType type() const;

    void setMin(float min);
    float min() const;

    void setMax(float max);
    float max() const;

    void setAutoAdjustRange(bool autoAdjust);
    bool isAutoAdjustRange() const;

    void setRange(float min, float max);

    void setLabelAutoRotation(float angle);
    float labelAutoRotation() const;

    void setTitleVisible(bool visible);
    bool isTitleVisible() const;

    void setTitleFixed(bool fixed);
    bool isTitleFixed() const;

Q_SIGNALS:
    void titleChanged(const QString &newTitle);
    void labelsChanged();
    void orientationChanged(QAbstract3DAxis::AxisOrientation orientation);
    void minChanged(float value);
    void maxChanged(float value);
    void rangeChanged(float min, float max);
    void autoAdjustRangeChanged(bool autoAdjust);
    Q_REVISION(1) void labelAutoRotationChanged(float angle);
    Q_REVISION(1) void titleVisibilityChanged(bool visible);
    Q_REVISION(1) void titleFixedChanged(bool fixed);

protected:
    QScopedPointer<QAbstract3DAxisPrivate> d_ptr;

private:
    Q_DISABLE_COPY(QAbstract3DAxis)

    friend class Abstract3DController;
    friend class Bars3DController;
    friend class QScatterDataProxyPrivate;
    friend class QSurfaceDataProxyPrivate;
};

QT_END_NAMESPACE_DATAVISUALIZATION

#endif
