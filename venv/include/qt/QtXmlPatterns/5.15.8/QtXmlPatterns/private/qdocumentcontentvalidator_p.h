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

#ifndef Patternist_DocumentContentValidator_H
#define Patternist_DocumentContentValidator_H

#include <private/qdynamiccontext_p.h>
#include <private/qexpression_p.h>
#include "qabstractxmlreceiver.h"

QT_BEGIN_NAMESPACE

namespace QPatternist
{
    /**
     * @short Receives QAbstractXmlReceiver events and validates that they are correct,
     * before sending them on to a second QAbstractXmlReceiver.
     *
     * @ingroup Patternist_xdm
     * @author Frans Englich <frans.englich@nokia.com>
     * @todo Escape data
     */
    class DocumentContentValidator : public QAbstractXmlReceiver
    {
    public:
        /**
         * DocumentContentValidator does not own @p receiver.
         */
        DocumentContentValidator(QAbstractXmlReceiver *const receiver,
                                 const DynamicContext::Ptr &context,
                                 const Expression::ConstPtr &expr);

        virtual void namespaceBinding(const QXmlName &nb);
        virtual void characters(const QStringRef &value);
        virtual void comment(const QString &value);

        virtual void startElement(const QXmlName &name);

        virtual void endElement();

        virtual void attribute(const QXmlName &name,
                               const QStringRef &value);

        virtual void processingInstruction(const QXmlName &name,
                                           const QString &value);

        virtual void item(const Item &item);

        virtual void startDocument();
        virtual void endDocument();
        virtual void atomicValue(const QVariant &value);
        virtual void startOfSequence();
        virtual void endOfSequence();

    private:
        QAbstractXmlReceiver *const      m_receiver;
        const DynamicContext::Ptr       m_context;
        const Expression::ConstPtr      m_expr;
        xsInteger                       m_elementDepth;
    };
}

QT_END_NAMESPACE

#endif
