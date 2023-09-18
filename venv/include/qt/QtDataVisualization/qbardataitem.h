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

#ifndef QBARDATAITEM_H
#define QBARDATAITEM_H

#include <QtDataVisualization/qdatavisualizationglobal.h>

QT_BEGIN_NAMESPACE_DATAVISUALIZATION

class QBarDataItemPrivate;

class QT_DATAVISUALIZATION_EXPORT QBarDataItem
{
public:
    QBarDataItem();
    QBarDataItem(float value);
    QBarDataItem(float value, float angle);
    QBarDataItem(const QBarDataItem &other);
    ~QBarDataItem();

    QBarDataItem &operator=(const QBarDataItem &other);

    inline void setValue(float val) { m_value = val; }
    inline float value() const { return m_value; }
    inline void setRotation(float angle) { m_angle = angle; }
    inline float rotation() const { return m_angle; }

protected:
    void createExtraData();

    QBarDataItemPrivate *d_ptr;

private:
    float m_value;
    float m_angle;
};

QT_END_NAMESPACE_DATAVISUALIZATION

#endif
