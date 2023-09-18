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

#ifndef CHARTELEMENT_H
#define CHARTELEMENT_H

#include <QtCharts/QChartGlobal>
#include <QtCore/QObject>
#include <QtCore/QRect>
#include <QGraphicsObject>
#include <QtCharts/private/qchartglobal_p.h>

QT_CHARTS_BEGIN_NAMESPACE

class ChartAnimator;
class ChartPresenter;
class ChartAnimation;
class ChartThemeManager;
class AbstractDomain;
class ChartDataSet;

class Q_CHARTS_PRIVATE_EXPORT ChartElement: public QGraphicsObject
{

public:
    explicit ChartElement(QGraphicsItem* item = 0);

    virtual ChartAnimation *animation() const { return 0; }
    virtual void setPresenter(ChartPresenter *presenter);
    ChartPresenter *presenter() const;
    virtual void setThemeManager(ChartThemeManager *manager);
    ChartThemeManager* themeManager() const;
    virtual void setDataSet(ChartDataSet *dataSet);
    ChartDataSet *dataSet() const;

private:
    ChartPresenter *m_presenter;
    ChartThemeManager *m_themeManager;
    ChartDataSet *m_dataSet;
};

QT_CHARTS_END_NAMESPACE

#endif
