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

#ifndef BOXWHISKERSDATA_P_H
#define BOXWHISKERSDATA_P_H

#include <QtCharts/QChartGlobal>
#include <QtCharts/private/qchartglobal_p.h>
#include <QtCore/QSizeF>

QT_CHARTS_BEGIN_NAMESPACE

class Q_CHARTS_PRIVATE_EXPORT BoxWhiskersData
{
public:
    BoxWhiskersData() :
        m_lowerExtreme(0.0),
        m_lowerQuartile(0.0),
        m_median(0.0),
        m_upperQuartile(0.0),
        m_upperExtreme(0.0),
        m_index(0),
        m_boxItems(0),
        m_maxX(0.0),
        m_minX(0.0),
        m_maxY(0.0),
        m_minY(0.0),
        m_seriesIndex(0),
        m_seriesCount(0)
    {
    }

    // Box related statistics
    qreal m_lowerExtreme;
    qreal m_lowerQuartile;
    qreal m_median;
    qreal m_upperQuartile;
    qreal m_upperExtreme;
    int   m_index;
    int   m_boxItems;

    // Domain boundaries, axis
    qreal m_maxX;
    qreal m_minX;
    qreal m_maxY;
    qreal m_minY;

    // Serieses related data
    int m_seriesIndex;
    int m_seriesCount;
};

QT_CHARTS_END_NAMESPACE

#endif // BOXWHISKERSDATA_P_H
