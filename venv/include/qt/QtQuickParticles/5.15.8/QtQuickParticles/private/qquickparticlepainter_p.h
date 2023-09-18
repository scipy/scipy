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

#ifndef PARTICLE_H
#define PARTICLE_H

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

#include <QObject>
#include <QDebug>
#include <QPair>
#include "qquickparticlesystem_p.h"

QT_BEGIN_NAMESPACE

class QQuickParticlePainter : public QQuickItem
{
    Q_OBJECT
    Q_PROPERTY(QQuickParticleSystem* system READ system WRITE setSystem NOTIFY systemChanged)
    Q_PROPERTY(QStringList groups READ groups WRITE setGroups NOTIFY groupsChanged)

    QML_NAMED_ELEMENT(ParticlePainter)
    QML_UNCREATABLE("Abstract type. Use one of the inheriting types instead.")

public: // data
    typedef QQuickParticleVarLengthArray<QQuickParticleGroupData::ID, 4> GroupIDs;

public:
    explicit QQuickParticlePainter(QQuickItem *parent = 0);
    //Data Interface to system
    void load(QQuickParticleData*);
    void reload(QQuickParticleData*);
    void setCount(int c);

    int count() const
    {
        return m_count;
    }

    void performPendingCommits();//Called from updatePaintNode
    QQuickParticleSystem* system() const
    {
        return m_system;
    }

    QStringList groups() const
    {
        return m_groups;
    }

    const GroupIDs &groupIds() const
    {
        if (m_groupIdsNeedRecalculation) {
            recalculateGroupIds();
        }
        return m_groupIds;
    }

    void itemChange(ItemChange, const ItemChangeData &) override;

Q_SIGNALS:
    void countChanged();
    void systemChanged(QQuickParticleSystem* arg);

    void groupsChanged(const QStringList &arg);

public Q_SLOTS:
    void setSystem(QQuickParticleSystem* arg);

    void setGroups(const QStringList &arg);

    void calcSystemOffset(bool resetPending = false);

private Q_SLOTS:
    virtual void sceneGraphInvalidated() {}

protected:
    /* Reset resets all your internal data structures. But anything attached to a particle should
       be in attached data. So reset + reloads should have no visible effect.
       ###Hunt down all cases where we do a complete reset for convenience and be more targeted
    */
    virtual void reset();

    void componentComplete() override;
    virtual void initialize(int gIdx, int pIdx){//Called from main thread
        Q_UNUSED(gIdx);
        Q_UNUSED(pIdx);
    }
    virtual void commit(int gIdx, int pIdx){//Called in Render Thread
        //###If you need to do something on size changed, check m_data size in this? Or we reset you every time?
        Q_UNUSED(gIdx);
        Q_UNUSED(pIdx);
    }

    QQuickParticleSystem* m_system;
    friend class QQuickParticleSystem;
    int m_count;
    bool m_pleaseReset;//Used by subclasses, but it's a nice optimization to know when stuff isn't going to matter.
    QPointF m_systemOffset;

    QQuickWindow *m_window;
    bool m_windowChanged;

private: // methods
    void recalculateGroupIds() const;

private: // data
    QStringList m_groups;
    QSet<QPair<int,int> > m_pendingCommits;
    mutable GroupIDs m_groupIds;
    mutable bool m_groupIdsNeedRecalculation;
};

QT_END_NAMESPACE
#endif // PARTICLE_H
