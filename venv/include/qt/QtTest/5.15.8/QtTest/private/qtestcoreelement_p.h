/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtTest module of the Qt Toolkit.
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

#ifndef QTESTCOREELEMENT_P_H
#define QTESTCOREELEMENT_P_H

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

#include <QtTest/private/qtestcorelist_p.h>
#include <QtTest/private/qtestelementattribute_p.h>

QT_BEGIN_NAMESPACE


template <class ElementType>
class QTestCoreElement: public QTestCoreList<ElementType>
{
    public:
        QTestCoreElement( int type = -1 );
        virtual ~QTestCoreElement();

        void addAttribute(const QTest::AttributeIndex index, const char *value);
        QTestElementAttribute *attributes() const;
        const char *attributeValue(QTest::AttributeIndex index) const;
        const char *attributeName(QTest::AttributeIndex index) const;
        const QTestElementAttribute *attribute(QTest::AttributeIndex index) const;

        const char *elementName() const;
        QTest::LogElementType elementType() const;

    private:
        QTestElementAttribute *listOfAttributes = nullptr;
        QTest::LogElementType type;
};

template<class ElementType>
QTestCoreElement<ElementType>::QTestCoreElement(int t)
    : type(QTest::LogElementType(t))
{
}

template<class ElementType>
QTestCoreElement<ElementType>::~QTestCoreElement()
{
    delete listOfAttributes;
}

template <class ElementType>
void QTestCoreElement<ElementType>::addAttribute(const QTest::AttributeIndex attributeIndex, const char *value)
{
    if (attributeIndex == -1 || attribute(attributeIndex))
        return;

    QTestElementAttribute *testAttribute = new QTestElementAttribute;
    testAttribute->setPair(attributeIndex, value);
    testAttribute->addToList(&listOfAttributes);
}

template <class ElementType>
QTestElementAttribute *QTestCoreElement<ElementType>::attributes() const
{
    return listOfAttributes;
}

template <class ElementType>
const char *QTestCoreElement<ElementType>::attributeValue(QTest::AttributeIndex index) const
{
    const QTestElementAttribute *attrb = attribute(index);
    if (attrb)
        return attrb->value();

    return nullptr;
}

template <class ElementType>
const char *QTestCoreElement<ElementType>::attributeName(QTest::AttributeIndex index) const
{
    const QTestElementAttribute *attrb = attribute(index);
    if (attrb)
        return attrb->name();

    return nullptr;
}

template <class ElementType>
const char *QTestCoreElement<ElementType>::elementName() const
{
    const char *xmlElementNames[] =
    {
        "property",
        "properties",
        "failure",
        "error",
        "testcase",
        "testsuite",
        "benchmark",
        "system-err"
    };

    if (type != QTest::LET_Undefined)
        return xmlElementNames[type];

    return nullptr;
}

template <class ElementType>
QTest::LogElementType QTestCoreElement<ElementType>::elementType() const
{
    return type;
}

template <class ElementType>
const QTestElementAttribute *QTestCoreElement<ElementType>::attribute(QTest::AttributeIndex index) const
{
    QTestElementAttribute *iterator = listOfAttributes;
    while (iterator) {
        if (iterator->index() == index)
            return iterator;

        iterator = iterator->nextElement();
    }

    return nullptr;
}

QT_END_NAMESPACE

#endif
