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

#ifndef QVALUE3DAXIS_H
#define QVALUE3DAXIS_H

#include <QtDataVisualization/qabstract3daxis.h>
#include <QtDataVisualization/qvalue3daxisformatter.h>

QT_BEGIN_NAMESPACE_DATAVISUALIZATION

class QValue3DAxisPrivate;

class QT_DATAVISUALIZATION_EXPORT QValue3DAxis : public QAbstract3DAxis
{
    Q_OBJECT
    Q_PROPERTY(int segmentCount READ segmentCount WRITE setSegmentCount NOTIFY segmentCountChanged)
    Q_PROPERTY(int subSegmentCount READ subSegmentCount WRITE setSubSegmentCount NOTIFY subSegmentCountChanged)
    Q_PROPERTY(QString labelFormat READ labelFormat WRITE setLabelFormat NOTIFY labelFormatChanged)
    Q_PROPERTY(QValue3DAxisFormatter* formatter READ formatter WRITE setFormatter NOTIFY formatterChanged REVISION 1)
    Q_PROPERTY(bool reversed READ reversed WRITE setReversed NOTIFY reversedChanged REVISION 1)

public:
    explicit QValue3DAxis(QObject *parent = nullptr);
    virtual ~QValue3DAxis();

    void setSegmentCount(int count);
    int segmentCount() const;

    void setSubSegmentCount(int count);
    int subSegmentCount() const;

    void setLabelFormat(const QString &format);
    QString labelFormat() const;

    void setFormatter(QValue3DAxisFormatter *formatter);
    QValue3DAxisFormatter *formatter() const;

    void setReversed(bool enable);
    bool reversed() const;

Q_SIGNALS:
    void segmentCountChanged(int count);
    void subSegmentCountChanged(int count);
    void labelFormatChanged(const QString &format);
    Q_REVISION(1) void formatterChanged(QValue3DAxisFormatter *formatter);
    Q_REVISION(1) void reversedChanged(bool enable);

protected:
    QValue3DAxisPrivate *dptr();
    const QValue3DAxisPrivate *dptrc() const;

private:
    Q_DISABLE_COPY(QValue3DAxis)
    friend class Abstract3DController;
    friend class Bars3DController;
    friend class Scatter3DController;
    friend class Surface3DController;
    friend class QValue3DAxisFormatterPrivate;
};

QT_END_NAMESPACE_DATAVISUALIZATION

#endif
