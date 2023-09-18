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

#ifndef Patternist_XsdReference_H
#define Patternist_XsdReference_H

#include <private/qxsdterm_p.h>

#include <QtXmlPatterns/QSourceLocation>

QT_BEGIN_NAMESPACE

namespace QPatternist
{
    /**
     * @short A helper class for element and group reference resolving.
     *
     * For easy resolving of element and group references, we have this class
     * that can be used as a place holder for the real element or group
     * object it is referring to.
     * So whenever the parser detects an element or group reference, it creates
     * a XsdReference and returns it instead of the XsdElement or XsdModelGroup.
     * During a later phase, the resolver will look for all XsdReferences
     * in the schema and will replace them with their referring XsdElement or
     * XsdModelGroup objects.
     *
     * @ingroup Patternist_schema
     * @author Tobias Koenig <tobias.koenig@nokia.com>
     */
    class XsdReference : public XsdTerm
    {
        public:
            typedef QExplicitlySharedDataPointer<XsdReference> Ptr;

            /**
             * Describes the type of the reference.
             */
            enum Type
            {
                Element,        ///< The reference points to an element.
                ModelGroup      ///< The reference points to a model group.
            };

            /**
             * Returns always @c true, used to avoid dynamic casts.
             */
            virtual bool isReference() const;

            /**
             * Sets the @p type of the reference.
             *
             * @see Type
             */
            void setType(Type type);

            /**
             * Returns the type of the reference.
             */
            Type type() const;

            /**
             * Sets the @p name of the referenced object.
             *
             * The name can either be a top-level element declaration
             * or a top-level group declaration.
             */
            void setReferenceName(const QXmlName &ame);

            /**
             * Returns the name of the referenced object.
             */
            QXmlName referenceName() const;

            /**
             * Sets the source @p location where the reference is located.
             */
            void setSourceLocation(const QSourceLocation &location);

            /**
             * Returns the source location where the reference is located.
             */
            QSourceLocation sourceLocation() const;

        private:
            Type            m_type;
            QXmlName        m_referenceName;
            QSourceLocation m_sourceLocation;
    };
}

QT_END_NAMESPACE

#endif
