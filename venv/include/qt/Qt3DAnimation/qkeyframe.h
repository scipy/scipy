/****************************************************************************
**
** Copyright (C) 2017 Klaralvdalens Datakonsult AB (KDAB).
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt3D module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:LGPL$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** GNU Lesser General Public License Usage
** Alternatively, this file may be used under the terms of the GNU Lesser
** General Public License version 3 as published by the Free Software
** Foundation and appearing in the file LICENSE.LGPL3 included in the
** packaging of this file. Please review the following information to
** ensure the GNU Lesser General Public License version 3 requirements
** will be met: https://www.gnu.org/licenses/lgpl-3.0.html.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 2.0 or (at your option) the GNU General
** Public license version 3 or any later version approved by the KDE Free
** Qt Foundation. The licenses are as published by the Free Software
** Foundation and appearing in the file LICENSE.GPL2 and LICENSE.GPL3
** included in the packaging of this file. Please review the following
** information to ensure the GNU General Public License requirements will
** be met: https://www.gnu.org/licenses/gpl-2.0.html and
** https://www.gnu.org/licenses/gpl-3.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef QT3DANIMATION_QKEYFRAME_H
#define QT3DANIMATION_QKEYFRAME_H

#include <QtGui/qvector2d.h>
#include <Qt3DAnimation/qt3danimation_global.h>

QT_BEGIN_NAMESPACE

namespace Qt3DAnimation {

class QKeyFrame
{
public:
    enum InterpolationType : quint8 {
        ConstantInterpolation,
        LinearInterpolation,
        BezierInterpolation
    };

    Q_DECL_CONSTEXPR QKeyFrame() Q_DECL_NOTHROW
        : m_coordinates()
        , m_leftControlPoint()
        , m_rightControlPoint()
        , m_interpolationType(BezierInterpolation)
    {
    }

    Q_DECL_CONSTEXPR explicit QKeyFrame(QVector2D coords) Q_DECL_NOTHROW
        : m_coordinates(coords)
        , m_leftControlPoint()
        , m_rightControlPoint()
        , m_interpolationType(LinearInterpolation)
    {
    }

    Q_DECL_CONSTEXPR explicit QKeyFrame(QVector2D coords,
                                        QVector2D lh,
                                        QVector2D rh) Q_DECL_NOTHROW
        : m_coordinates(coords)
        , m_leftControlPoint(lh)
        , m_rightControlPoint(rh)
        , m_interpolationType(BezierInterpolation)
    {
    }

    void setCoordinates(QVector2D coords) Q_DECL_NOTHROW
    {
        m_coordinates = coords;
    }

    Q_DECL_CONSTEXPR QVector2D coordinates() const Q_DECL_NOTHROW
    {
        return m_coordinates;
    }

    void setLeftControlPoint(QVector2D lh) Q_DECL_NOTHROW
    {
        m_leftControlPoint = lh;
    }

    Q_DECL_CONSTEXPR QVector2D leftControlPoint() const Q_DECL_NOTHROW
    {
        return m_leftControlPoint;
    }

    void setRightControlPoint(QVector2D rh) Q_DECL_NOTHROW
    {
        m_rightControlPoint = rh;
    }

    Q_DECL_CONSTEXPR QVector2D rightControlPoint() const Q_DECL_NOTHROW
    {
        return m_rightControlPoint;
    }

    void setInterpolationType(InterpolationType interp) Q_DECL_NOTHROW
    {
        m_interpolationType = interp;
    }

    Q_DECL_CONSTEXPR InterpolationType interpolationType() const Q_DECL_NOTHROW
    {
        return m_interpolationType;
    }

    friend inline bool operator==(const QKeyFrame &, const QKeyFrame &) Q_DECL_NOTHROW;
    friend inline bool operator!=(const QKeyFrame &, const QKeyFrame &) Q_DECL_NOTHROW;

private:
    QVector2D m_coordinates;
    QVector2D m_leftControlPoint;
    QVector2D m_rightControlPoint;
    InterpolationType m_interpolationType;
};

QT3D_DECLARE_TYPEINFO(Qt3DAnimation, QKeyFrame, Q_PRIMITIVE_TYPE)

inline bool operator==(const QKeyFrame &lhs, const QKeyFrame &rhs) Q_DECL_NOTHROW
{
    if (lhs.m_interpolationType != rhs.m_interpolationType)
        return false;

    if (lhs.m_interpolationType == QKeyFrame::BezierInterpolation) {
        return lhs.m_coordinates == rhs.m_coordinates &&
               lhs.m_leftControlPoint == rhs.m_leftControlPoint &&
               lhs.m_rightControlPoint == rhs.m_rightControlPoint;
    }

    return lhs.m_coordinates == rhs.m_coordinates;
}

inline bool operator!=(const QKeyFrame &lhs, const QKeyFrame &rhs) Q_DECL_NOTHROW
{
    return !(lhs == rhs);
}

} // namespace Qt3DAnimation

QT_END_NAMESPACE

#endif // QT3DANIMATION_QKEYFRAME_H
