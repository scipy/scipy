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

#ifndef QXMLFORMATTER_H
#define QXMLFORMATTER_H

#include <QtXmlPatterns/QXmlSerializer>

QT_BEGIN_NAMESPACE


class QIODevice;
class QTextCodec;
class QXmlQuery;
class QXmlFormatterPrivate;

class Q_XMLPATTERNS_EXPORT QXmlFormatter : public QXmlSerializer
{
public:
    QXmlFormatter(const QXmlQuery &query,
                  QIODevice *outputDevice);

    void characters(const QStringRef &value) override;
    void comment(const QString &value) override;
    void startElement(const QXmlName &name) override;
    void endElement() override;

    void attribute(const QXmlName &name,
                   const QStringRef &value) override;
    void processingInstruction(const QXmlName &name,
                               const QString &value) override;
    void atomicValue(const QVariant &value) override;
    void startDocument() override;
    void endDocument() override;
    void startOfSequence() override;
    void endOfSequence() override;

    int indentationDepth() const;
    void setIndentationDepth(int depth);

    /* The members below are internal, not part of the public API, and
     * unsupported. Using them leads to undefined behavior. */
    void item(const QPatternist::Item &item) override;
private:
    inline void startFormattingContent();
    Q_DECLARE_PRIVATE(QXmlFormatter)
};

QT_END_NAMESPACE

#endif
