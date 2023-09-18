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

#ifndef Patternist_XsdStateMachineBuilder_H
#define Patternist_XsdStateMachineBuilder_H

#include <private/qxsdparticle_p.h>
#include <private/qxsdstatemachine_p.h>
#include <private/qxsdterm_p.h>

#include <QtCore/QExplicitlySharedDataPointer>
#include <QtCore/QList>

QT_BEGIN_NAMESPACE

namespace QPatternist
{
    /**
     * @short A helper class to build up validation state machines.
     *
     * @ingroup Patternist_schema
     * @author Tobias Koenig <tobias.koenig@nokia.com>
     */
    class XsdStateMachineBuilder : public QSharedData
    {
        public:
            typedef QExplicitlySharedDataPointer<XsdStateMachineBuilder> Ptr;

            enum Mode
            {
                CheckingMode,
                ValidatingMode
            };

            /**
             * Creates a new state machine builder.
             *
             * @param machine The state machine it should work on.
             * @param namePool The name pool used by all schema components.
             * @param mode The mode the machine shall be build for.
             */
            XsdStateMachineBuilder(XsdStateMachine<XsdTerm::Ptr> *machine, const NamePool::Ptr &namePool, Mode mode = CheckingMode);

            /**
             * Resets the state machine.
             *
             * @returns The initial end state.
             */
            XsdStateMachine<XsdTerm::Ptr>::StateId reset();

            /**
             * Prepends a start state to the given @p state.
             * That is needed to allow the conversion of the state machine from a FSA to a DFA.
             */
            XsdStateMachine<XsdTerm::Ptr>::StateId addStartState(XsdStateMachine<XsdTerm::Ptr>::StateId state);

            /**
             * Creates the state machine for the given @p particle that should have the
             * given @p endState.
             *
             * @returns The new start state.
             */
            XsdStateMachine<XsdTerm::Ptr>::StateId buildParticle(const XsdParticle::Ptr &particle, XsdStateMachine<XsdTerm::Ptr>::StateId endState);

            /**
             * Creates the state machine for the given @p term that should have the
             * given @p endState.
             *
             * @returns The new start state.
             */
            XsdStateMachine<XsdTerm::Ptr>::StateId buildTerm(const XsdTerm::Ptr &term, XsdStateMachine<XsdTerm::Ptr>::StateId endState);

            /**
             * Returns a hash that maps each term that appears inside @p particle, to the particle it belongs.
             *
             * @note These information are used by XsdParticleChecker to check particle inheritance.
             */
            static QHash<XsdTerm::Ptr, XsdParticle::Ptr> particleLookupMap(const XsdParticle::Ptr &particle);

        private:
            XsdStateMachine<XsdTerm::Ptr> *m_stateMachine;
            NamePool::Ptr                  m_namePool;
            Mode                           m_mode;
    };
}

QT_END_NAMESPACE

#endif
