/****************************************************************************
**
** Copyright (C) 2015 The Qt Company Ltd.
** Contact: http://www.qt.io/licensing/
**
** This file is part of the QtScript module of the Qt Toolkit.
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

#ifndef QSCRIPTPROGRAM_H
#define QSCRIPTPROGRAM_H

#include <QtCore/qsharedpointer.h>

#include <QtCore/qstring.h>
#include <QtScript/qtscriptglobal.h>

QT_BEGIN_NAMESPACE


class QScriptProgramPrivate;
class Q_SCRIPT_EXPORT QScriptProgram
{
public:
    QScriptProgram();
    QScriptProgram(const QString &sourceCode,
                   const QString fileName = QString(),
                   int firstLineNumber = 1);
    QScriptProgram(const QScriptProgram &other);
    ~QScriptProgram();

    QScriptProgram &operator=(const QScriptProgram &other);

    bool isNull() const;

    QString sourceCode() const;
    QString fileName() const;
    int firstLineNumber() const;

    bool operator==(const QScriptProgram &other) const;
    bool operator!=(const QScriptProgram &other) const;

private:
    QExplicitlySharedDataPointer<QScriptProgramPrivate> d_ptr;
    Q_DECLARE_PRIVATE(QScriptProgram)
};

QT_END_NAMESPACE

#endif // QSCRIPTPROGRAM_H
