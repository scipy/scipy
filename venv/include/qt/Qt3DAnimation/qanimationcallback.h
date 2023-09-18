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

#ifndef QT3DANIMATION_QANIMATIONCALLBACK_H
#define QT3DANIMATION_QANIMATIONCALLBACK_H

#include <Qt3DAnimation/qt3danimation_global.h>
#include <Qt3DCore/qnode.h>

QT_BEGIN_NAMESPACE

namespace Qt3DAnimation {

class Q_3DANIMATIONSHARED_EXPORT QAnimationCallback
{
public:
    enum Flag {
        OnOwningThread = 0x0,
        OnThreadPool = 0x01
    };
    Q_DECLARE_FLAGS(Flags, Flag)

    virtual ~QAnimationCallback() { }

    virtual void valueChanged(const QVariant &value) = 0;
};

} // namespace Qt3DAnimation

Q_DECLARE_OPERATORS_FOR_FLAGS(Qt3DAnimation::QAnimationCallback::Flags)

QT_END_NAMESPACE

#endif // QT3DANIMATION_QANIMATIONCALLBACK_H
