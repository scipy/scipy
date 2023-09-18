/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt Charts module of the Qt Toolkit.
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

//  W A R N I N G
//  -------------
//
// This file is not part of the Qt Chart API.  It exists purely as an
// implementation detail.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.

#ifndef DECLARATIVE_MARGINS_H
#define DECLARATIVE_MARGINS_H

#include <QtCharts/QChartGlobal>
#include <QtCore/QObject>
#include <QtCore/QMargins>
#include <private/declarativechartglobal_p.h>

QT_CHARTS_BEGIN_NAMESPACE

class Q_QMLCHARTS_PRIVATE_EXPORT DeclarativeMargins : public QObject, public QMargins
{
    Q_OBJECT
    Q_PROPERTY(int top READ top WRITE setTop NOTIFY topChanged)
    Q_PROPERTY(int bottom READ bottom WRITE setBottom NOTIFY bottomChanged)
    Q_PROPERTY(int left READ left WRITE setLeft NOTIFY leftChanged)
    Q_PROPERTY(int right READ right WRITE setRight NOTIFY rightChanged)

public:
    explicit DeclarativeMargins(QObject *parent = 0);
    void setTop(int top);
    void setBottom(int bottom);
    void setLeft(int left);
    void setRight(int right);

Q_SIGNALS:
    void topChanged(int top, int bottom, int left, int right);
    void bottomChanged(int top, int bottom, int left, int right);
    void leftChanged(int top, int bottom, int left, int right);
    void rightChanged(int top, int bottom, int left, int right);
};

QT_CHARTS_END_NAMESPACE

#endif // DECLARATIVE_MARGINS_H
