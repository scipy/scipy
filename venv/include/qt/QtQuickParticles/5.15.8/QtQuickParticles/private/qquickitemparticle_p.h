/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtQuick module of the Qt Toolkit.
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

#ifndef ITEMPARTICLE_H
#define ITEMPARTICLE_H

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
#include "qquickparticlepainter_p.h"
#include <QPointer>
#include <QSet>
#include <private/qquickanimation_p_p.h>
QT_BEGIN_NAMESPACE

class QQuickItemParticleAttached;

class QQuickItemParticle : public QQuickParticlePainter
{
    Q_OBJECT
    Q_PROPERTY(bool fade READ fade WRITE setFade NOTIFY fadeChanged)
    Q_PROPERTY(QQmlComponent* delegate READ delegate WRITE setDelegate NOTIFY delegateChanged)
    QML_NAMED_ELEMENT(ItemParticle)
    QML_ATTACHED(QQuickItemParticleAttached)
public:
    explicit QQuickItemParticle(QQuickItem *parent = 0);
    ~QQuickItemParticle();

    bool fade() const { return m_fade; }

    QSGNode *updatePaintNode(QSGNode *, UpdatePaintNodeData *) override;

    static QQuickItemParticleAttached *qmlAttachedProperties(QObject *object);
    QQmlComponent* delegate() const
    {
        return m_delegate;
    }

Q_SIGNALS:
    void fadeChanged();

    void delegateChanged(QQmlComponent* arg);

public Q_SLOTS:
    //TODO: Add a follow mode, where moving the delegate causes the logical particle to go with it?
    void freeze(QQuickItem* item);
    void unfreeze(QQuickItem* item);
    void take(QQuickItem* item,bool prioritize=false);//take by modelparticle
    void give(QQuickItem* item);//give from modelparticle

    void setFade(bool arg){if (arg == m_fade) return; m_fade = arg; Q_EMIT fadeChanged();}
    void setDelegate(QQmlComponent* arg)
    {
        if (m_delegate != arg) {
            m_delegate = arg;
            Q_EMIT delegateChanged(arg);
        }
    }

protected:
    void reset() override;
    void commit(int gIdx, int pIdx) override;
    void initialize(int gIdx, int pIdx) override;
    void prepareNextFrame();
private:
    void processDeletables();
    void tick(int time = 0);
    QSet<QQuickItem* > m_deletables;
    QList<QQuickItem* > m_managed;
    bool m_fade;

    QList<QQuickItem*> m_pendingItems;
    QList<int> m_available;
    QSet<QQuickItem*> m_stasis;
    qreal m_lastT;
    int m_activeCount;
    QQmlComponent* m_delegate;

    typedef QTickAnimationProxy<QQuickItemParticle, &QQuickItemParticle::tick> Clock;
    Clock *clock;
};

class QQuickItemParticleAttached : public QObject
{
    Q_OBJECT
    Q_PROPERTY(QQuickItemParticle* particle READ particle CONSTANT);
public:
    QQuickItemParticleAttached(QObject* parent)
        : QObject(parent), m_mp(0), m_parentItem(nullptr)
    {;}
    QQuickItemParticle* particle() const { return m_mp; }
    void detach(){Q_EMIT detached();}
    void attach(){Q_EMIT attached();}
private:
    QQuickItemParticle* m_mp;
    QPointer<QQuickItem> m_parentItem;
    friend class QQuickItemParticle;
Q_SIGNALS:
    void detached();
    void attached();
};

QT_END_NAMESPACE

#endif // ITEMPARTICLE_H
