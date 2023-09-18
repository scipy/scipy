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

#ifndef Patternist_XsdDocumentation_H
#define Patternist_XsdDocumentation_H

#include <private/qanytype_p.h>
#include <private/qanyuri_p.h>
#include <private/qderivedstring_p.h>
#include <private/qnamedschemacomponent_p.h>

QT_BEGIN_NAMESPACE

namespace QPatternist
{
    /**
     * @short Represents a XSD documentation object.
     *
     * This class represents the <em>documentation</em> component of an <em>annotation</em> object
     * of a XML schema as described <a href="http://www.w3.org/TR/xmlschema11-1/#cAnnotations">here</a>.
     *
     * @ingroup Patternist_schema
     * @author Tobias Koenig <tobias.koenig@nokia.com>
     */
    class XsdDocumentation : public NamedSchemaComponent
    {
        public:
            typedef QExplicitlySharedDataPointer<XsdDocumentation> Ptr;
            typedef QList<XsdDocumentation::Ptr> List;

            /**
             * Creates a new documentation object.
             */
            XsdDocumentation();

            /**
             * Destroys the documentation object.
             */
            ~XsdDocumentation();

            /**
             * Sets the @p source of the documentation.
             *
             * The source points to an URL that contains more
             * information.
             */
            void setSource(const AnyURI::Ptr &source);

            /**
             * Returns the source of the documentation.
             */
            AnyURI::Ptr source() const;

            /**
             * Sets the @p language of the documentation.
             */
            void setLanguage(const DerivedString<TypeLanguage>::Ptr &language);

            /**
             * Returns the language of the documentation.
             */
            DerivedString<TypeLanguage>::Ptr language() const;

            /**
             * Sets the @p content of the documentation.
             *
             * The content can be of abritrary type.
             */
            void setContent(const QString &content);

            /**
             * Returns the content of the documentation.
             */
            QString content() const;

        private:
            AnyURI::Ptr                      m_source;
            DerivedString<TypeLanguage>::Ptr m_language;
            QString                          m_content;
    };
}

QT_END_NAMESPACE

#endif
