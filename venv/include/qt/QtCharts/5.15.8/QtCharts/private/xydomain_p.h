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

#ifndef XYDOMAIN_H
#define XYDOMAIN_H
#include <private/abstractdomain_p.h>
#include <QtCharts/private/qchartglobal_p.h>
#include <QtCore/QRectF>
#include <QtCore/QSizeF>

QT_CHARTS_BEGIN_NAMESPACE

class Q_CHARTS_PRIVATE_EXPORT XYDomain: public AbstractDomain
{
    Q_OBJECT
public:
    explicit XYDomain(QObject *object = 0);
    virtual ~XYDomain();

    DomainType type(){ return AbstractDomain::XYDomain;}

    void setRange(qreal minX, qreal maxX, qreal minY, qreal maxY);

    friend bool Q_AUTOTEST_EXPORT operator== (const XYDomain &Domain1, const XYDomain &Domain2);
    friend bool Q_AUTOTEST_EXPORT operator!= (const XYDomain &Domain1, const XYDomain &Domain2);
    friend QDebug Q_AUTOTEST_EXPORT operator<<(QDebug dbg, const XYDomain &AbstractDomain);

    void zoomIn(const QRectF &rect);
    void zoomOut(const QRectF &rect);
    void move(qreal dx, qreal dy);

    QPointF calculateGeometryPoint(const QPointF &point, bool &ok) const;
    QPointF calculateDomainPoint(const QPointF &point) const;
    QVector<QPointF> calculateGeometryPoints(const QVector<QPointF> &vector) const;
};

QT_CHARTS_END_NAMESPACE

#endif // XYDOMAIN_H
