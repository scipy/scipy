/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtScxml module of the Qt Toolkit.
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

#ifndef QSCXMLSTATEMACHINEINFO_H
#define QSCXMLSTATEMACHINEINFO_H

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

#include <QtScxml/qscxmlglobals.h>
#include <QtCore/qobject.h>

QT_BEGIN_NAMESPACE

class QScxmlStateMachine;
class QScxmlStateMachineInfoPrivate;

class Q_SCXML_EXPORT QScxmlStateMachineInfo: public QObject
{
    Q_OBJECT

public: // types
    typedef int StateId;
    typedef int TransitionId;

    static const StateId InvalidStateId = -1;
    static const TransitionId InvalidTransitionId = -1;

    enum StateType : int {
        InvalidState = -1,
        NormalState = 0,
        ParallelState = 1,
        FinalState = 2,
        ShallowHistoryState = 3,
        DeepHistoryState = 4
    };

    enum TransitionType : int {
        InvalidTransition = -1,
        InternalTransition = 0,
        ExternalTransition = 1,
        SyntheticTransition = 2
    };

public: // methods
    QScxmlStateMachineInfo(QScxmlStateMachine *stateMachine);

    QScxmlStateMachine *stateMachine() const;

    QVector<StateId> allStates() const;
    QVector<TransitionId> allTransitions() const;
    QString stateName(int stateId) const;
    StateId stateParent(StateId stateId) const;
    StateType stateType(int stateId) const;
    QVector<StateId> stateChildren(StateId stateId) const;
    TransitionId initialTransition(StateId stateId) const;
    TransitionType transitionType(TransitionId transitionId) const;
    StateId transitionSource(TransitionId transitionId) const;
    QVector<StateId> transitionTargets(TransitionId transitionId) const;
    QVector<QString> transitionEvents(TransitionId transitionId) const;
    QVector<StateId> configuration() const;

Q_SIGNALS:
    void statesEntered(const QVector<QScxmlStateMachineInfo::StateId> &states);
    void statesExited(const QVector<QScxmlStateMachineInfo::StateId> &states);
    void transitionsTriggered(const QVector<QScxmlStateMachineInfo::TransitionId> &transitions);

private:
    Q_DECLARE_PRIVATE(QScxmlStateMachineInfo)
};

QT_END_NAMESPACE

#endif // QSCXMLSTATEMACHINEINFO_H
