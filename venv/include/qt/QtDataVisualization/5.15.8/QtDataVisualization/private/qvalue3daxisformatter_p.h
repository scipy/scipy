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

//
//  W A R N I N G
//  -------------
//
// This file is not part of the QtDataVisualization API.  It exists purely as an
// implementation detail.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.

#ifndef QVALUE3DAXISFORMATTER_P_H
#define QVALUE3DAXISFORMATTER_P_H

#include "datavisualizationglobal_p.h"
#include "qvalue3daxisformatter.h"
#include "utils_p.h"
#include <QtCore/QLocale>

QT_BEGIN_NAMESPACE_DATAVISUALIZATION

class QValue3DAxis;

class QValue3DAxisFormatterPrivate : public QObject
{
    Q_OBJECT

public:
    QValue3DAxisFormatterPrivate(QValue3DAxisFormatter *q);
    virtual ~QValue3DAxisFormatterPrivate();

    void recalculate();
    void doRecalculate();
    void populateCopy(QValue3DAxisFormatter &copy);
    void doPopulateCopy(QValue3DAxisFormatterPrivate &copy);

    QString stringForValue(qreal value, const QString &format);
    float positionAt(float value) const;
    float valueAt(float position) const;

    void setAxis(QValue3DAxis *axis);
    void markDirty(bool labelsChange);

public Q_SLOTS:
    void markDirtyNoLabelChange();

protected:
    QValue3DAxisFormatter *q_ptr;

    bool m_needsRecalculate;

    float m_min;
    float m_max;
    float m_rangeNormalizer;

    QVector<float> m_gridPositions;
    QVector<float> m_subGridPositions;
    QVector<float> m_labelPositions;
    QStringList m_labelStrings;

    QValue3DAxis *m_axis;

    QString m_previousLabelFormat;
    QByteArray m_labelFormatArray;
    Utils::ParamType m_preparsedParamType;

    bool m_allowNegatives;
    bool m_allowZero;

    QLocale m_locale;
    QString m_formatPreStr;
    QString m_formatPostStr;
    int m_formatPrecision;
    char m_formatSpec;
    bool m_cLocaleInUse;

    friend class QValue3DAxisFormatter;
};

QT_END_NAMESPACE_DATAVISUALIZATION

#endif
