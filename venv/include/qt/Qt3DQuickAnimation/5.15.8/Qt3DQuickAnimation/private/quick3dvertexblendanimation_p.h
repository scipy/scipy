/****************************************************************************
**
** Copyright (C) 2017 The Qt Company Ltd.
** Contact: http://www.qt.io/licensing/
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

#ifndef QT3DANIMATION_QUICK_QUICK3DVERTEXBLENDANIMATION_P_H
#define QT3DANIMATION_QUICK_QUICK3DVERTEXBLENDANIMATION_P_H

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

#include <QtQml/qqmllist.h>

#include <Qt3DAnimation/qvertexblendanimation.h>

#include <private/qt3dquickanimation_global_p.h>

QT_BEGIN_NAMESPACE

namespace Qt3DAnimation {
namespace Quick {

class Q_3DQUICKANIMATIONSHARED_PRIVATE_EXPORT QQuick3DVertexBlendAnimation : public QObject
{
    Q_OBJECT
    Q_PROPERTY(QQmlListProperty<Qt3DAnimation::QMorphTarget> morphTargets READ morphTargets)

public:
    QQuick3DVertexBlendAnimation(QObject *parent = nullptr);

    inline QVertexBlendAnimation *parentVertexBlendAnimation() const
    {
        return qobject_cast<QVertexBlendAnimation *>(parent());
    }

    QQmlListProperty<Qt3DAnimation::QMorphTarget> morphTargets();

private:

    static void appendMorphTarget(QQmlListProperty<Qt3DAnimation::QMorphTarget> *list,
                                  Qt3DAnimation::QMorphTarget *morphTarget);
    static Qt3DAnimation::QMorphTarget *morphTargetAt(
                                  QQmlListProperty<Qt3DAnimation::QMorphTarget> *list,
                                  int index);
    static int morphTargetCount(QQmlListProperty<Qt3DAnimation::QMorphTarget> *list);
    static void clearMorphTargets(QQmlListProperty<Qt3DAnimation::QMorphTarget> *list);
};

} // namespace Quick
} // namespace Qt3DAnimation

QT_END_NAMESPACE

#endif // QT3DANIMATION_QUICK_QUICK3DMORPHINGANIMATION_P_H
