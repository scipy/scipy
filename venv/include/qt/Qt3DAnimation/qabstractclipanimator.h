/****************************************************************************
**
** Copyright (C) 2016 Klaralvdalens Datakonsult AB (KDAB).
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

#ifndef QT3DANIMATION_QABSTRACTCLIPANIMATOR_H
#define QT3DANIMATION_QABSTRACTCLIPANIMATOR_H

#include <Qt3DAnimation/qt3danimation_global.h>
#include <Qt3DCore/qcomponent.h>

QT_BEGIN_NAMESPACE

namespace Qt3DAnimation {

class QAnimationClip;
class QChannelMapper;
class QClock;
class QAbstractClipAnimatorPrivate;

class Q_3DANIMATIONSHARED_EXPORT QAbstractClipAnimator : public Qt3DCore::QComponent
{
    Q_OBJECT
    Q_PROPERTY(bool running READ isRunning WRITE setRunning NOTIFY runningChanged)
    Q_PROPERTY(int loops READ loopCount WRITE setLoopCount NOTIFY loopCountChanged)
    Q_PROPERTY(Qt3DAnimation::QChannelMapper *channelMapper READ channelMapper WRITE setChannelMapper NOTIFY channelMapperChanged)
    Q_PROPERTY(Qt3DAnimation::QClock *clock READ clock WRITE setClock NOTIFY clockChanged)
    Q_PROPERTY(float normalizedTime READ normalizedTime WRITE setNormalizedTime NOTIFY normalizedTimeChanged)

public:
    enum Loops { Infinite = -1 };
    Q_ENUM(Loops)

    ~QAbstractClipAnimator();

    bool isRunning() const;
    Qt3DAnimation::QChannelMapper *channelMapper() const;
    int loopCount() const;
    Qt3DAnimation::QClock *clock() const;
    float normalizedTime() const;

public Q_SLOTS:
    void setRunning(bool running);
    void setChannelMapper(Qt3DAnimation::QChannelMapper *channelMapper);
    void setLoopCount(int loops);
    void setClock(Qt3DAnimation::QClock *clock);
    void setNormalizedTime(float timeFraction);

    void start();
    void stop();

Q_SIGNALS:
    void runningChanged(bool running);
    void channelMapperChanged(Qt3DAnimation::QChannelMapper *channelMapper);
    void loopCountChanged(int loops);
    void clockChanged(Qt3DAnimation::QClock *clock);
    void normalizedTimeChanged(float index);

protected:
    explicit QAbstractClipAnimator(Qt3DCore::QNode *parent = nullptr);
    QAbstractClipAnimator(QAbstractClipAnimatorPrivate &dd, Qt3DCore::QNode *parent = nullptr);

    // TODO Unused remove in Qt6
    void sceneChangeEvent(const Qt3DCore::QSceneChangePtr &change) override;
private:
    Q_DECLARE_PRIVATE(QAbstractClipAnimator)
};

} // namespace Qt3DAnimation

QT_END_NAMESPACE

#endif // QT3DANIMATION_QABSTRACTCLIPANIMATOR_H
