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

#ifndef QLOGVALUE3DAXISFORMATTER_H
#define QLOGVALUE3DAXISFORMATTER_H

#include <QtDataVisualization/qvalue3daxisformatter.h>

QT_BEGIN_NAMESPACE_DATAVISUALIZATION

class QLogValue3DAxisFormatterPrivate;

class QT_DATAVISUALIZATION_EXPORT QLogValue3DAxisFormatter : public QValue3DAxisFormatter
{
    Q_OBJECT

    Q_PROPERTY(qreal base READ base WRITE setBase NOTIFY baseChanged)
    Q_PROPERTY(bool autoSubGrid READ autoSubGrid WRITE setAutoSubGrid NOTIFY autoSubGridChanged)
    Q_PROPERTY(bool showEdgeLabels READ showEdgeLabels WRITE setShowEdgeLabels NOTIFY showEdgeLabelsChanged)

protected:
    explicit QLogValue3DAxisFormatter(QLogValue3DAxisFormatterPrivate *d, QObject *parent = nullptr);
public:
    explicit QLogValue3DAxisFormatter(QObject *parent = nullptr);
    virtual ~QLogValue3DAxisFormatter();

    void setBase(qreal base);
    qreal base() const;
    void setAutoSubGrid(bool enabled);
    bool autoSubGrid() const;
    void setShowEdgeLabels(bool enabled);
    bool showEdgeLabels() const;

Q_SIGNALS:
    void baseChanged(qreal base);
    void autoSubGridChanged(bool enabled);
    void showEdgeLabelsChanged(bool enabled);

protected:
    virtual QValue3DAxisFormatter *createNewInstance() const;
    virtual void recalculate();
    virtual float positionAt(float value) const;
    virtual float valueAt(float position) const;
    virtual void populateCopy(QValue3DAxisFormatter &copy) const;

    QLogValue3DAxisFormatterPrivate *dptr();
    const QLogValue3DAxisFormatterPrivate *dptrc() const;

private:
    Q_DISABLE_COPY(QLogValue3DAxisFormatter)

    friend class QLogValue3DAxisFormatterPrivate;
};

QT_END_NAMESPACE_DATAVISUALIZATION

#endif
