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

#ifndef Patternist_TokenSource_H
#define Patternist_TokenSource_H

#include <private/qfunctionargument_p.h>
#include <private/qitem_p.h>
#include <private/qitemtype_p.h>
#include <private/qtokenvalue_p.h>
#include <private/qparsercontext_p.h>
#include <private/qquerytransformparser_p.h>

QT_BEGIN_NAMESPACE

template<typename T> class QQueue;

namespace QPatternist
{
    /**
     * @short Base class for components that needs to return tokens.
     *
     * TokenSource represents a stream of Token instances. The end
     * is reached when readNext() returns a Token constructed with
     * END_OF_FILE.
     *
     * @see <a href="http://www.w3.org/TR/xquery-xpath-parsing/">Building a
     * Tokenizer for XPath or XQuery</a>
     * @author Frans Englich <frans.englich@nokia.com>
     */
    class TokenSource : public QSharedData
    {
    public:
        /**
         * typedef for the enum Bison generates that contains
         * the token symbols.
         */
        typedef XPathtokentype TokenType;

        /**
         * Represents a token by carrying its name and value.
         */
        class Token
        {
        public:
            /**
             * Constructs an invalid Token. This default constructor
             * is need in Qt's container classes.
             */
            inline Token() {}
            inline Token(const TokenType t) : type(t) {}
            inline Token(const TokenType t, const QString &val) : type(t), value(val) {}

            bool hasError() const
            {
                return type == T_ERROR;
            }

            TokenType type;
            QString value;
        };

        typedef QExplicitlySharedDataPointer<TokenSource> Ptr;
        typedef QQueue<Ptr> Queue;

        /**
         * The C++ compiler cannot synthesize it when we use the
         * Q_DISABLE_COPY() macro.
         */
        inline TokenSource()
        {
        }

        virtual ~TokenSource();

        /**
         * @returns the next token.
         */
        virtual Token nextToken(XPATHLTYPE *const sourceLocator) = 0;

    private:
        Q_DISABLE_COPY(TokenSource)
    };
}

QT_END_NAMESPACE

#endif
