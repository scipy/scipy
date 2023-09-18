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

#ifndef QSSGPERFTIMER_H
#define QSSGPERFTIMER_H

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

#include <QtCore/QVector>
#include <QtCore/qmutex.h>
#include <QtCore/qhash.h>
#include <QtCore/qelapsedtimer.h>

QT_BEGIN_NAMESPACE

struct QSSGTimerEntry;

class Q_QUICK3DUTILS_EXPORT QSSGPerfTimer
{
    Q_DISABLE_COPY(QSSGPerfTimer)

public:
    struct Key {
        const char *id;
    };

    struct Entry {
        Entry(const QString &id)
            : tag(id)
        {}
        Entry() = default;

        void update(qint64 elapsed);
        void reset();
        QString toString(quint32 nFrames) const;

        quint32 count = 0;
        qint64 totalTime = 0;
        qint64 maxTime = 0;
        QString tag;
    };

private:
    bool m_isEnabled = false;
    int frameCount = 0;
    QMutex mutex;
    QHash<Key, Entry> entries;

public:
    QSSGPerfTimer();
    ~QSSGPerfTimer();

    // amount is in counter frequency units
    void update(const char *inTag, qint64 inAmount);

    // Dump current summation of timer data.
    void dump();
    void reset();

    int newFrame() { return ++frameCount; }

    void setEnabled(bool b) { m_isEnabled = b; }
    bool isEnabled() const { return m_isEnabled; }
};

struct QSSGStackPerfTimer
{
    QSSGPerfTimer *m_timer;
    QElapsedTimer elapsedTimer;
    const char *m_id;

    QSSGStackPerfTimer(QSSGPerfTimer *timer, const char *inId)
        : m_timer(timer), m_id(inId)
    {
        if (timer->isEnabled())
            elapsedTimer.start();
    }

    ~QSSGStackPerfTimer()
    {
        if (m_timer->isEnabled()) {
            qint64 elapsed = elapsedTimer.nsecsElapsed();
            m_timer->update(m_id, elapsed);
        }
    }
};

QT_END_NAMESPACE

#endif // QSSGPERFTIMER_H
