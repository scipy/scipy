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

#ifndef CHARTTHEMEBLUECERULEAN_P_H
#define CHARTTHEMEBLUECERULEAN_P_H

#include <private/charttheme_p.h>
#include <QtCharts/private/qchartglobal_p.h>

QT_CHARTS_BEGIN_NAMESPACE

class Q_CHARTS_PRIVATE_EXPORT ChartThemeBlueCerulean: public ChartTheme
{
public:
    ChartThemeBlueCerulean() : ChartTheme(QChart::ChartThemeBlueCerulean) {
        // Series colors
        m_seriesColors << QRgb(0xc7e85b);
        m_seriesColors << QRgb(0x1cb54f);
        m_seriesColors << QRgb(0x5cbf9b);
        m_seriesColors << QRgb(0x009fbf);
        m_seriesColors << QRgb(0xee7392);
        m_seriesGradients = ChartThemeManager::generateSeriesGradients(m_seriesColors);

        // Background
        QLinearGradient backgroundGradient(0.5, 0.0, 0.5, 1.0);
        backgroundGradient.setColorAt(0.0, QRgb(0x056189));
        backgroundGradient.setColorAt(1.0, QRgb(0x101a31));
        backgroundGradient.setCoordinateMode(QGradient::ObjectBoundingMode);
        m_chartBackgroundGradient = backgroundGradient;

        // Axes and other
        m_labelBrush = QBrush(QRgb(0xffffff));
        m_axisLinePen = QPen(QRgb(0xd6d6d6));
        m_axisLinePen.setWidth(2);
        m_gridLinePen = QPen(QRgb(0x84a2b0));
        m_gridLinePen.setWidth(1);
        m_minorGridLinePen = QPen(QRgb(0x84a2b0));
        m_minorGridLinePen.setWidth(1);
        m_minorGridLinePen.setStyle(Qt::DashLine);
        m_backgroundShades = BackgroundShadesNone;
        m_outlinePen = QPen(QRgb(0xebebeb));
        m_outlinePen.setWidthF(2.0);
    }
};

QT_CHARTS_END_NAMESPACE

#endif
