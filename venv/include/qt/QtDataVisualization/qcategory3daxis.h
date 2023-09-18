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

#ifndef QCATEGORY3DAXIS_H
#define QCATEGORY3DAXIS_H

#include <QtDataVisualization/qabstract3daxis.h>

QT_BEGIN_NAMESPACE_DATAVISUALIZATION

class QCategory3DAxisPrivate;

class QT_DATAVISUALIZATION_EXPORT QCategory3DAxis : public QAbstract3DAxis
{
    Q_OBJECT
    Q_PROPERTY(QStringList labels READ labels WRITE setLabels NOTIFY labelsChanged)

public:
    explicit QCategory3DAxis(QObject *parent = nullptr);
    virtual ~QCategory3DAxis();

    void setLabels(const QStringList &labels);
    QStringList labels() const;

Q_SIGNALS:
    void labelsChanged();

protected:
    QCategory3DAxisPrivate *dptr();

private:
    Q_DISABLE_COPY(QCategory3DAxis)
    friend class Bars3DController;
};

QT_END_NAMESPACE_DATAVISUALIZATION

#endif
