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

#ifndef Patternist_XsdTerm_H
#define Patternist_XsdTerm_H

#include <private/qnamedschemacomponent_p.h>
#include <private/qxsdannotated_p.h>

QT_BEGIN_NAMESPACE

namespace QPatternist
{
    /**
     * @short A base class for all particles of a model group.
     *
     * This class is the base class for all particles of a model group
     * as the <em>element</em>, <em>group</em> or <em>any</em> tag, it is not supposed to
     * be instantiated directly.
     *
     * @see <a href="http://www.w3.org/Submission/2004/SUBM-xmlschema-api-20040309/xml-schema-api.html#Interface-XSTerm">XML Schema API reference</a>
     * @ingroup Patternist_schema
     * @author Tobias Koenig <tobias.koenig@nokia.com>
     */
    class XsdTerm : public NamedSchemaComponent, public XsdAnnotated
    {
        public:
            typedef QExplicitlySharedDataPointer<XsdTerm> Ptr;

            /**
             * Returns @c true if the term is an element, @c false otherwise.
             */
            virtual bool isElement() const;

            /**
             * Returns @c true if the term is a model group (group tag), @c false otherwise.
             */
            virtual bool isModelGroup() const;

            /**
             * Returns @c true if the term is a wildcard (any tag), @c false otherwise.
             */
            virtual bool isWildcard() const;

            /**
             * Returns @c true if the term is a reference, @c false otherwise.
             *
             * @note The reference term is only used internally as helper during type resolving.
             */
            virtual bool isReference() const;

        protected:
            /**
             * This constructor only exists to ensure this class is subclassed.
             */
            inline XsdTerm() {};
    };
}

QT_END_NAMESPACE

#endif
