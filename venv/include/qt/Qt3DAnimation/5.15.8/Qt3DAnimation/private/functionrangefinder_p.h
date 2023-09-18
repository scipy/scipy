/****************************************************************************
**
** Copyright (C) 2017 Klaralvdalens Datakonsult AB (KDAB).
** Contact: http://www.qt-project.org/legal
**
** This file is part of the Qt3D module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:LGPL3$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see http://www.qt.io/terms-conditions. For further
** information use the contact form at http://www.qt.io/contact-us.
**
** GNU Lesser General Public License Usage
** Alternatively, this file may be used under the terms of the GNU Lesser
** General Public License version 3 as published by the Free Software
** Foundation and appearing in the file LICENSE.LGPLv3 included in the
** packaging of this file. Please review the following information to
** ensure the GNU Lesser General Public License version 3 requirements
** will be met: https://www.gnu.org/licenses/lgpl.html.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 2.0 or later as published by the Free
** Software Foundation and appearing in the file LICENSE.GPL included in
** the packaging of this file. Please review the following information to
** ensure the GNU General Public License version 2.0 requirements will be
** met: http://www.gnu.org/licenses/gpl-2.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef QT3DANIMATION_ANIMATION_FUNCTIONRANGEFINDER_P_H
#define QT3DANIMATION_ANIMATION_FUNCTIONRANGEFINDER_P_H

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API.  It exists for the convenience
// of other Qt classes.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.
//

#include <QtCore/qvector.h>

#include <cmath>
#include <cstdlib>

QT_BEGIN_NAMESPACE

namespace Qt3DAnimation {
namespace Animation {

class Q_AUTOTEST_EXPORT FunctionRangeFinder
{
public:
    FunctionRangeFinder(const QVector<float> &x);

    inline int findLowerBound(float x) const { return m_correlated ? hunt(x) : locate(x); }

    int rangeSize() const { return m_rangeSize; }
    void setRangeSize(int rangeSize) { m_rangeSize = rangeSize; }

    bool isAscending() const { return m_ascending; }
    void setAscending(bool ascending) { m_ascending = ascending; }

    int correlationThreshold() const { return m_correlationThreshold; }
    void updateAutomaticCorrelationThreshold()
    {
        m_correlationThreshold = std::max(1, int(std::pow(float(m_x.size()), 0.25)));
    }

private:
    int locate(float x) const;
    int hunt(float x) const;

    const QVector<float> &m_x;
    mutable int m_previousLowerBound;
    mutable bool m_correlated;
    int m_rangeSize;
    int m_correlationThreshold;
    bool m_ascending;
};

} // namespace Animation
} // namespace Qt3DAnimation

QT_END_NAMESPACE

#endif // QT3DANIMATION_ANIMATION_FUNCTIONRANGEFINDER_P_H
