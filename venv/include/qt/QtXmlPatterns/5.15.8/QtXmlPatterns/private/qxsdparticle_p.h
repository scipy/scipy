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

#ifndef Patternist_XsdParticle_H
#define Patternist_XsdParticle_H

#include <private/qnamedschemacomponent_p.h>
#include <private/qxsdterm_p.h>

#include <QtCore/QList>

QT_BEGIN_NAMESPACE

namespace QPatternist
{
    /**
     * @short Represents a XSD particle object.
     *
     * This class represents the <em>particle</em> object of a XML schema as described
     * <a href="http://www.w3.org/TR/xmlschema11-1/#cParticles">here</a>.
     *
     * It contains information about the number of occurrence and a reference to
     * either an <em>element</em> object, a <em>group</em> object or an <em>any</em> object.
     *
     * @see <a href="http://www.w3.org/Submission/2004/SUBM-xmlschema-api-20040309/xml-schema-api.html#Interface-XSParticle">XML Schema API reference</a>
     * @ingroup Patternist_schema
     * @author Tobias Koenig <tobias.koenig@nokia.com>
     */
    class XsdParticle : public NamedSchemaComponent
    {
        public:
            typedef QExplicitlySharedDataPointer<XsdParticle> Ptr;
            typedef QList<XsdParticle::Ptr> List;

            /**
             * Creates a new particle object.
             */
            XsdParticle();

            /**
             * Sets the minimum @p occurrence of the particle.
             *
             * @see <a href="http://www.w3.org/TR/xmlschema11-1/#p-min_occurs">Minimum Occurrence Definition</a>
             */
            void setMinimumOccurs(unsigned int occurrence);

            /**
             * Returns the minimum occurrence of the particle.
             */
            unsigned int minimumOccurs() const;

            /**
             * Sets the maximum @p occurrence of the particle.
             *
             * @see <a href="http://www.w3.org/TR/xmlschema11-1/#p-max_occurs">Maximum Occurrence Definition</a>
             */
            void setMaximumOccurs(unsigned int occurrence);

            /**
             * Returns the maximum occurrence of the particle.
             *
             * @note This value has only a meaning if maximumOccursUnbounded is @c false.
             */
            unsigned int maximumOccurs() const;

            /**
             * Sets whether the maximum occurrence of the particle is unbounded.
             *
             * @see <a href="http://www.w3.org/TR/xmlschema11-1/#p-max_occurs">Maximum Occurrence Definition</a>
             */
            void setMaximumOccursUnbounded(bool unbounded);

            /**
             * Returns whether the maximum occurrence of the particle is unbounded.
             */
            bool maximumOccursUnbounded() const;

            /**
             * Sets the @p term of the particle.
             *
             * The term can be an element, a model group or an element wildcard.
             *
             * @see <a href="http://www.w3.org/TR/xmlschema11-1/#p-term">Term Definition</a>
             */
            void setTerm(const XsdTerm::Ptr &term);

            /**
             * Returns the term of the particle.
             */
            XsdTerm::Ptr term() const;

        private:
            unsigned int m_minimumOccurs;
            unsigned int m_maximumOccurs;
            bool         m_maximumOccursUnbounded;
            XsdTerm::Ptr m_term;
    };
}

QT_END_NAMESPACE

#endif
