/****************************************************************************
**
** Copyright (C) 2017 The Qt Company Ltd.
** Contact: http://www.qt.io/licensing/
**
** This file is part of the Qt Quick Controls 2 module of the Qt Toolkit.
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

#ifndef QQUICKANIMATEDNODE_P_H
#define QQUICKANIMATEDNODE_P_H

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

#include <QtQuick/qsgnode.h>
#include <QtCore/qelapsedtimer.h>
#include <QtQuickControls2/private/qtquickcontrols2global_p.h>

QT_BEGIN_NAMESPACE

class QQuickItem;
class QQuickWindow;

class Q_QUICKCONTROLS2_PRIVATE_EXPORT QQuickAnimatedNode : public QObject, public QSGTransformNode
{
    Q_OBJECT

public:
    explicit QQuickAnimatedNode(QQuickItem *target);

    bool isRunning() const;

    int currentTime() const;
    void setCurrentTime(int time);

    int duration() const;
    void setDuration(int duration);

    enum LoopCount { Infinite = -1 };

    int loopCount() const;
    void setLoopCount(int count);

    virtual void sync(QQuickItem *target);

    QQuickWindow *window() const;

    // must be called from sync() or updatePaintNode()
    void start(int duration = 0);
    void restart();
    void stop();

Q_SIGNALS:
    void started();
    void stopped();

protected:
    virtual void updateCurrentTime(int time);

private Q_SLOTS:
    void advance();
    void update();

private:
    bool m_running = false;
    int m_duration = 0;
    int m_loopCount = 1;
    int m_currentTime = 0;
    int m_currentLoop = 0;
    QElapsedTimer m_timer;
    QQuickWindow *m_window = nullptr;
};

QT_END_NAMESPACE

#endif // QQUICKANIMATEDNODE_P_H
