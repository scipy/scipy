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

#ifndef QABSTRACTXMLPULLPROVIDER_H
#define QABSTRACTXMLPULLPROVIDER_H

#include <QtCore/QtGlobal>

QT_BEGIN_NAMESPACE

class QXmlItem;
class QXmlName;
class QString;
class QVariant;
template<typename Key, typename Value> class QHash;

namespace QPatternist
{
    class AbstractXmlPullProviderPrivate;

    class AbstractXmlPullProvider
    {
    public:
        AbstractXmlPullProvider();
        virtual ~AbstractXmlPullProvider();

        enum Event
        {
            StartOfInput            = 1,
            AtomicValue             = 1 << 1,
            StartDocument           = 1 << 2,
            EndDocument             = 1 << 3,
            StartElement            = 1 << 4,
            EndElement              = 1 << 5,
            Text                    = 1 << 6,
            ProcessingInstruction   = 1 << 7,
            Comment                 = 1 << 8,
            Attribute               = 1 << 9,
            Namespace               = 1 << 10,
            EndOfInput              = 1 << 11
        };

        virtual Event next() = 0;
        virtual Event current() const = 0;
        virtual QXmlName name() const = 0;
        virtual QVariant atomicValue() const = 0;
        virtual QString stringValue() const = 0;

        virtual QHash<QXmlName, QString> attributes() = 0;
        virtual QHash<QXmlName, QXmlItem> attributeItems() = 0;

        /* *** The functions below are internal. */
    private:
        Q_DISABLE_COPY(AbstractXmlPullProvider)
    };
}

QT_END_NAMESPACE

#endif
