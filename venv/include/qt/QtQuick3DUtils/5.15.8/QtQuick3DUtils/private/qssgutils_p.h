/****************************************************************************
**
** Copyright (C) 2008-2012 NVIDIA Corporation.
** Copyright (C) 2019 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of Qt Quick 3D.
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

#ifndef QSSGUTILS_H
#define QSSGUTILS_H

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API.  It exists purely as an
// implementation detail.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.
//

#include <QtQuick3DUtils/private/qtquick3dutilsglobal_p.h>
#include <QtQuick3DUtils/private/qssgdataref_p.h>

#include <QtGui/QVector2D>
#include <QtGui/QVector3D>
#include <QtGui/QQuaternion>
#include <QtGui/QMatrix3x3>
#include <QtGui/QMatrix4x4>

#include <QtCore/qdebug.h>
#include <QtCore/QString>
#include <QtCore/qloggingcategory.h>
#include <QtCore/QIODevice>
#include <QtCore/qmath.h>

QT_BEGIN_NAMESPACE

namespace aux {
Q_DECL_CONSTEXPR inline float translateBrightness(float brightness) { return brightness * .01f; }
Q_DECL_CONSTEXPR inline float translateConstantAttenuation(float attenuation) { return attenuation; }
template<int MINATTENUATION = 0, int MAXATTENUATION = 1000>
Q_DECL_CONSTEXPR inline float translateLinearAttenuation(float attenuation) { return qBound(float(MINATTENUATION), attenuation, float(MAXATTENUATION)) * .01f; }
template<int MINATTENUATION = 0, int MAXATTENUATION = 1000>
Q_DECL_CONSTEXPR inline float translateQuadraticAttenuation(float attenuation) { return qBound(float(MINATTENUATION), attenuation, float(MAXATTENUATION)) * .0001f; }
}

namespace vec2 {
float Q_QUICK3DUTILS_EXPORT magnitude(const QVector2D &v);
}

namespace vec3 {
inline QVector3D minimum(const QVector3D &v1, const QVector3D &v2) Q_DECL_NOTHROW { return { qMin(v1.x(), v2.x()), qMin(v1.y(), v2.y()), qMin(v1.z(), v2.z()) }; }
inline QVector3D maximum(const QVector3D &v1, const QVector3D &v2) Q_DECL_NOTHROW { return { qMax(v1.x(), v2.x()), qMax(v1.y(), v2.y()), qMax(v1.z(), v2.z()) }; }
bool Q_QUICK3DUTILS_EXPORT isFinite(const QVector3D &v);
float Q_QUICK3DUTILS_EXPORT magnitude(const QVector3D &v);
float Q_QUICK3DUTILS_EXPORT magnitudeSquared(const QVector3D &v);
float Q_QUICK3DUTILS_EXPORT normalize(QVector3D &v);
}

namespace mat33 {
QVector3D Q_QUICK3DUTILS_EXPORT transform(const QMatrix3x3 &m, const QVector3D &v);
QMatrix3x3 Q_QUICK3DUTILS_EXPORT getInverse(const QMatrix3x3 &m);
}

namespace mat44 {
QMatrix3x3 Q_QUICK3DUTILS_EXPORT getUpper3x3(const QMatrix4x4 &m);
void Q_QUICK3DUTILS_EXPORT normalize(QMatrix4x4 &m);
QVector3D Q_QUICK3DUTILS_EXPORT rotate(const QMatrix4x4 &m, const QVector3D &v);
QVector4D Q_QUICK3DUTILS_EXPORT rotate(const QMatrix4x4 &m, const QVector4D &v);
QVector3D Q_QUICK3DUTILS_EXPORT transform(const QMatrix4x4 &m, const QVector3D &v);
QVector4D Q_QUICK3DUTILS_EXPORT transform(const QMatrix4x4 &m, const QVector4D &v);
QVector3D Q_QUICK3DUTILS_EXPORT getPosition(const QMatrix4x4 &m);
QVector3D Q_QUICK3DUTILS_EXPORT getScale(const QMatrix4x4 &m);

inline void flip(QMatrix4x4 &matrix)
{
    // Flip between left-handed and right-handed orientation
    float *writePtr(matrix.data());
    // rotation conversion
    writePtr[0 * 4 + 2] *= -1;
    writePtr[1 * 4 + 2] *= -1;
    writePtr[2 * 4 + 0] *= -1;
    writePtr[2 * 4 + 1] *= -1;
    // translation conversion
    writePtr[3 * 4 + 2] *= -1;
}

}

namespace quant {
bool Q_QUICK3DUTILS_EXPORT isFinite(const QQuaternion &q);

float Q_QUICK3DUTILS_EXPORT magnitude(const QQuaternion &q);

bool Q_QUICK3DUTILS_EXPORT isSane(const QQuaternion &q);

bool Q_QUICK3DUTILS_EXPORT isUnit(const QQuaternion &q);

QVector3D Q_QUICK3DUTILS_EXPORT rotated(const QQuaternion &q, const QVector3D &v);

QVector3D Q_QUICK3DUTILS_EXPORT inverseRotated(const QQuaternion &q, const QVector3D &v);
}

template<typename TDataType>
QSSGDataRef<TDataType> PtrAtOffset(quint8 *baseData, quint32 offset, quint32 byteSize)
{
    return QSSGDataRef<TDataType>(byteSize ? reinterpret_cast<TDataType *>(baseData + offset) : nullptr,
                                    byteSize / sizeof(TDataType));
}

Q_QUICK3DUTILS_EXPORT const char *nonNull(const char *src);

inline QVector3D degToRad(const QVector3D &v) {
    return QVector3D(qDegreesToRadians(v.x()), qDegreesToRadians(v.y()), qDegreesToRadians(v.z()));
}

inline QVector3D radToDeg(const QVector3D &v) {
    return QVector3D(qRadiansToDegrees(v.x()), qRadiansToDegrees(v.y()), qRadiansToDegrees(v.z()));
}

QT_END_NAMESPACE

#endif // QSSGUTILS_H
