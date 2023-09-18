/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtXmlPatterns module of the Qt Toolkit.
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

/*
 * NOTE: This file is included by qxsdstatemachine_p.h
 * if you need some includes, put them in qxsdstatemachine_p.h (outside of the namespace)
 */

template <typename TransitionType>
XsdStateMachine<TransitionType>::XsdStateMachine()
    : m_counter(50),
      m_lastTransition()
{
}

template <typename TransitionType>
XsdStateMachine<TransitionType>::XsdStateMachine(const NamePool::Ptr &namePool)
    : m_namePool(namePool),
      m_counter(50),
      m_lastTransition()
{
}

template <typename TransitionType>
typename XsdStateMachine<TransitionType>::StateId XsdStateMachine<TransitionType>::addState(StateType type)
{
#ifndef QT_NO_DEBUG
    // make sure we don't have two start states
    if (type == StartState) {
        for (auto it = m_states.cbegin(), end = m_states.cend(); it != end; ++it) {
            Q_ASSERT(it.value() != StartState && it.value() != StartEndState);
        }
    }
#endif // QT_NO_DEBUG

    // reserve new state id
    const StateId id = ++m_counter;
    m_states.insert(id, type);

    // if it is a start state, we make it to our current state
    if (type == StartState || type == StartEndState)
        m_currentState = id;

    return id;
}

template <typename TransitionType>
void XsdStateMachine<TransitionType>::addTransition(StateId start, TransitionType transition, StateId end)
{
    QHash<TransitionType, QVector<StateId> > &hash = m_transitions[start];
    QVector<StateId> &states = hash[transition];
    if (!states.contains(end))
        states.append(end);
}

template <typename TransitionType>
void XsdStateMachine<TransitionType>::addEpsilonTransition(StateId start, StateId end)
{
    QVector<StateId> &states = m_epsilonTransitions[start];
    states.append(end);
}

template <typename TransitionType>
void XsdStateMachine<TransitionType>::reset()
{
    // reset the machine to the start state
    auto it = m_states.cbegin();
    auto end = m_states.cend();
    for ( ; it != end; ++it) {
        if (it.value() == StartState || it.value() == StartEndState) {
            m_currentState = it.key();
            return;
        }
    }

    Q_ASSERT(false);
}

template <typename TransitionType>
void XsdStateMachine<TransitionType>::clear()
{
    m_states.clear();
    m_transitions.clear();
    m_epsilonTransitions.clear();
    m_currentState = -1;
    m_counter = 50;
}

template <typename TransitionType>
bool XsdStateMachine<TransitionType>::proceed(TransitionType transition)
{
    // check that we are not in an invalid state
    if (!m_transitions.contains(m_currentState)) {
        return false;
    }

    // fetch the transition entry for the current state
    const QHash<TransitionType, QVector<StateId> > &entry = m_transitions[m_currentState];
    if (entry.contains(transition)) { // is there an transition for the given input?
        m_currentState = entry.value(transition).first();
        m_lastTransition = transition;
        return true;
    } else {
        return false;
    }
}

template <typename TransitionType>
QList<TransitionType> XsdStateMachine<TransitionType>::possibleTransitions() const
{
    // check that we are not in an invalid state
    if (!m_transitions.contains(m_currentState)) {
        return QList<TransitionType>();
    }

    // fetch the transition entry for the current state
    const QHash<TransitionType, QVector<StateId> > &entry = m_transitions[m_currentState];

    return entry.keys();
}

template <typename TransitionType>
template <typename InputType>
bool XsdStateMachine<TransitionType>::proceed(InputType input)
{
    // check that we are not in an invalid state
    if (!m_transitions.contains(m_currentState)) {
        return false;
    }

    // fetch the transition entry for the current state
    const QHash<TransitionType, QVector<StateId> > &entry = m_transitions[m_currentState];
    auto it = entry.cbegin();
    auto end = entry.cend();
    for ( ; it != end; ++it) {
        if (inputEqualsTransition(input, it.key())) {
            m_currentState = it.value().first();
            m_lastTransition = it.key();
            return true;
        }
    }

    return false;
}

template <typename TransitionType>
template <typename InputType>
bool XsdStateMachine<TransitionType>::inputEqualsTransition(InputType input, TransitionType transition) const
{
    Q_UNUSED(input);
    Q_UNUSED(transition);

    return false;
}

template <typename TransitionType>
bool XsdStateMachine<TransitionType>::inEndState() const
{
    // check if current state is an end state
    return (m_states.value(m_currentState) == StartEndState || m_states.value(m_currentState) == EndState);
}

template <typename TransitionType>
TransitionType XsdStateMachine<TransitionType>::lastTransition() const
{
    return m_lastTransition;
}

template <typename TransitionType>
typename XsdStateMachine<TransitionType>::StateId XsdStateMachine<TransitionType>::startState() const
{
    auto it = m_states.cbegin();
    auto end = m_states.cend();
    for ( ; it != end; ++it) {
        if (it.value() == StartState || it.value() == StartEndState)
            return it.key();
    }

    Q_ASSERT(false); // should never be reached
    return -1;
}

template <typename TransitionType>
QString XsdStateMachine<TransitionType>::transitionTypeToString(TransitionType type) const
{
    Q_UNUSED(type)

    return QString();
}

template <typename TransitionType>
bool XsdStateMachine<TransitionType>::outputGraph(QIODevice *device, const QString &graphName) const
{
    if (!device->isOpen()) {
        qWarning("device must be open for writing");
        return false;
    }

    QByteArray graph;
    QTextStream s(&graph);

    s << "digraph " << graphName << " {\n";
    s << "  mindist = 2.0\n";

    // draw edges
    for (auto it = m_transitions.cbegin(), end = m_transitions.cend(); it != end; ++it) {

        for (auto it2 = it.value().cbegin(), end = it.value().cend(); it2 != end; ++it2) {
            for (int i = 0; i < it2.value().count(); ++i)
                s << "  " << it.key() << " -> " << it2.value().at(i) << " [label=\"" << transitionTypeToString(it2.key()) << "\"]\n";
        }
    }

    for (auto it = m_epsilonTransitions.cbegin(), end = m_epsilonTransitions.cend(); it != end; ++it) {
        const QVector<StateId> states = it.value();
        for (int i = 0; i < states.count(); ++i)
            s << "  " << it.key() << " -> " << states.at(i) << " [label=\"&#949;\"]\n";
    }

    // draw node info
    for (auto it = m_states.cbegin(), end = m_states.cend(); it != end; ++it) {

        QString style;
        if (it.value() == StartState) {
            style = QLatin1String("shape=circle, style=filled, color=blue");
        } else if (it.value() == StartEndState) {
            style = QLatin1String("shape=doublecircle, style=filled, color=blue");
        } else if (it.value() == InternalState) {
            style = QLatin1String("shape=circle, style=filled, color=red");
        } else if (it.value() == EndState) {
            style = QLatin1String("shape=doublecircle, style=filled, color=green");
        }

        s << "  " << it.key() << " [" << style << "]\n";
    }

    s << "}\n";

    s.flush();

    if (device->write(graph) == -1)
        return false;

    return true;
}


template <typename TransitionType>
typename XsdStateMachine<TransitionType>::StateId XsdStateMachine<TransitionType>::dfaStateForNfaState(QSet<StateId> nfaState,
                                                                                                       QList< QPair<QSet<StateId>, StateId> > &stateTable,
                                                                                                       XsdStateMachine<TransitionType> &dfa) const
{
    // check whether we have the given state in our lookup table
    // already, in that case simply return it
    for (int i = 0; i < stateTable.count(); ++i) {
        if (stateTable.at(i).first == nfaState)
            return stateTable.at(i).second;
    }

    // check if the NFA state set contains a Start or End
    // state, in that case our new DFA state will be a
    // Start or End state as well
    StateType type = InternalState;
    bool hasStartState = false;
    bool hasEndState = false;
    for (const StateId state : qAsConst(nfaState)) {
        if (m_states.value(state) == EndState) {
            hasEndState = true;
        } else if (m_states.value(state) == StartState) {
            hasStartState = true;
        }
    }
    if (hasStartState) {
        if (hasEndState)
            type = StartEndState;
        else
            type = StartState;
    } else if (hasEndState) {
        type = EndState;
    }

    // create the new DFA state
    const StateId dfaState = dfa.addState(type);

    // add the new DFA state to the lookup table
    stateTable.append(qMakePair<QSet<StateId>, StateId>(nfaState, dfaState));

    return dfaState;
}

template <typename TransitionType>
XsdStateMachine<TransitionType> XsdStateMachine<TransitionType>::toDFA() const
{
    XsdStateMachine<TransitionType> dfa(m_namePool);
    dfa.m_counter = 100;
    QList< QPair< QSet<StateId>, StateId> > table;
    QList< QSet<StateId> > isMarked;

    // search the start state as the algorithm starts with it...
    StateId startState = -1;
    auto it = m_states.cbegin();
    auto end = m_states.cend();
    for ( ; it != end; ++it) {
        if (it.value() == StartState) {
            startState = it.key();
            break;
        }
    }
    Q_ASSERT(startState != -1);

    // our list of state set that still have to be processed
    QList< QSet<StateId> > workStates;

    // add the start state to the list of to processed state sets
    auto firstDfaState = epsilonClosure(QSet<StateId>() << startState);
    workStates.append(firstDfaState);
    isMarked.append(firstDfaState);

    while (!workStates.isEmpty()) { // as long as there are state sets to process left
        // enqueue set of states
        const QSet<StateId> states = workStates.takeFirst();

        // select a list of all inputs that are possible for
        // the 'states' set
        QList<TransitionType> input;

        for (const StateId state : states)
            input << m_transitions.value(state).keys();

        // get the state in DFA that corresponds to the 'states' set in the NFA
        const StateId dfaBegin = dfaStateForNfaState(states, table, dfa);

        for (int i = 0; i < input.count(); ++i) { // for each possible input
            // retrieve the states that can  be reached from the 'states' set by the
            // given input or by epsilon transition
            const QSet<StateId> followStates = epsilonClosure(move(states, input.at(i)));

            // get the state in DFA that corresponds to the 'followStates' set in the NFA
            const StateId dfaEnd = dfaStateForNfaState(followStates, table, dfa);

            // adds a new transition to the DFA that corresponds to the transitions between
            // 'states' and 'followStates' in the NFA
            dfa.addTransition(dfaBegin, input.at(i), dfaEnd);

            // add the 'followStates' to the list of to be processed state sets
            if (!isMarked.contains(followStates)) {
                workStates.append(followStates);
                isMarked.append(followStates); // only needs to be processed once
            }
        }
    }

    return dfa;
}

template <typename TransitionType>
QHash<typename XsdStateMachine<TransitionType>::StateId, typename XsdStateMachine<TransitionType>::StateType> XsdStateMachine<TransitionType>::states() const
{
    return m_states;
}

template <typename TransitionType>
QHash<typename XsdStateMachine<TransitionType>::StateId, QHash<TransitionType, QVector<typename XsdStateMachine<TransitionType>::StateId> > > XsdStateMachine<TransitionType>::transitions() const
{
    return m_transitions;
}
