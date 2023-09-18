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

#ifndef CANDLESTICKDATA_P_H
#define CANDLESTICKDATA_P_H

#include <QtCharts/QCandlestickSeries>
#include <QtCharts/private/qchartglobal_p.h>

QT_CHARTS_BEGIN_NAMESPACE

class Q_CHARTS_PRIVATE_EXPORT CandlestickData
{
public:
    CandlestickData() :
        m_timestamp(0.0),
        m_open(0.0),
        m_high(0.0),
        m_low(0.0),
        m_close(0.0),
        m_index(0),
        m_maxX(0.0),
        m_minX(0.0),
        m_maxY(0.0),
        m_minY(0.0),
        m_series(nullptr),
        m_seriesIndex(0),
        m_seriesCount(0)
    {
    }

    // Candlestick related statistics
    qreal m_timestamp;
    qreal m_open;
    qreal m_high;
    qreal m_low;
    qreal m_close;
    int m_index;

    // Domain boundaries
    qreal m_maxX;
    qreal m_minX;
    qreal m_maxY;
    qreal m_minY;

    // Series related data
    QCandlestickSeries *m_series;
    int m_seriesIndex;
    int m_seriesCount;
};

QT_CHARTS_END_NAMESPACE

#endif // CANDLESTICKDATA_P_H
