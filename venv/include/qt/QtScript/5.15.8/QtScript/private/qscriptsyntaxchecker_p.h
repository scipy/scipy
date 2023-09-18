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

#ifndef QSCRIPTSYNTAXCHECKER_H
#define QSCRIPTSYNTAXCHECKER_H

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

#include <QtCore/qstring.h>

#include "qscriptgrammar_p.h"

#include <stdlib.h>

QT_BEGIN_NAMESPACE

namespace QScript {

class Lexer;

class SyntaxChecker: protected QScriptGrammar
{
public:
    enum State {
        Error,
        Intermediate,
        Valid,
    };

    struct Result {
        Result(State s, int ln, int col, const QString &msg)
            : state(s), errorLineNumber(ln), errorColumnNumber(col),
              errorMessage(msg) {}
        State state;
        int errorLineNumber;
        int errorColumnNumber;
        QString errorMessage;
    };

    SyntaxChecker();
    ~SyntaxChecker();

    Result checkSyntax(const QString &code);

protected:
    bool automatic(QScript::Lexer *lexer, int token) const;
    inline void reallocateStack();

protected:
    int tos;
    int stack_size;
    int *state_stack;
};

inline void SyntaxChecker::reallocateStack()
{
    if (! stack_size)
        stack_size = 128;
    else
        stack_size <<= 1;

    state_stack = reinterpret_cast<int*> (realloc(state_stack, stack_size * sizeof(int)));
}

} // namespace QScript

QT_END_NAMESPACE

#endif
