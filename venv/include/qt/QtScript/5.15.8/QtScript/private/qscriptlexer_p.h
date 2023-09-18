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

#ifndef QSCRIPTLEXER_P_H
#define QSCRIPTLEXER_P_H

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

QT_BEGIN_NAMESPACE

class QScriptEnginePrivate;
class QScriptNameIdImpl;

namespace QScript {

class Lexer
{
public:
    Lexer(QScriptEnginePrivate *eng);
    ~Lexer();

    void setCode(const QString &c, int lineno);
    int lex();

    int currentLineNo() const { return yylineno; }
    int currentColumnNo() const { return yycolumn; }

    int startLineNo() const { return startlineno; }
    int startColumnNo() const { return startcolumn; }

    int endLineNo() const { return currentLineNo(); }
    int endColumnNo() const
    { int col = currentColumnNo(); return (col > 0) ? col - 1 : col; }

    bool prevTerminator() const { return terminator; }

    enum State { Start,
                 Identifier,
                 InIdentifier,
                 InSingleLineComment,
                 InMultiLineComment,
                 InNum,
                 InNum0,
                 InHex,
                 InOctal,
                 InDecimal,
                 InExponentIndicator,
                 InExponent,
                 Hex,
                 Octal,
                 Number,
                 String,
                 Eof,
                 InString,
                 InEscapeSequence,
                 InHexEscape,
                 InUnicodeEscape,
                 Other,
                 Bad };

    enum Error {
        NoError,
        IllegalCharacter,
        UnclosedStringLiteral,
        IllegalEscapeSequence,
        IllegalUnicodeEscapeSequence,
        UnclosedComment,
        IllegalExponentIndicator,
        IllegalIdentifier
    };

    enum ParenthesesState {
        IgnoreParentheses,
        CountParentheses,
        BalancedParentheses
    };

    enum RegExpBodyPrefix {
        NoPrefix,
        EqualPrefix
    };

    bool scanRegExp(RegExpBodyPrefix prefix = NoPrefix);

    QScriptNameIdImpl *pattern;
    int flags;

    State lexerState() const
        { return state; }

    QString errorMessage() const
        { return errmsg; }
    void setErrorMessage(const QString &err)
        { errmsg = err; }
    void setErrorMessage(const char *err)
        { setErrorMessage(QString::fromLatin1(err)); }

    Error error() const
        { return err; }
    void clearError()
        { err = NoError; }

private:
    QScriptEnginePrivate *driver;
    int yylineno;
    bool done;
    char *buffer8;
    QChar *buffer16;
    uint size8, size16;
    uint pos8, pos16;
    bool terminator;
    bool restrKeyword;
    // encountered delimiter like "'" and "}" on last run
    bool delimited;
    int stackToken;

    State state;
    void setDone(State s);
    uint pos;
    void shift(uint p);
    int lookupKeyword(const char *);

    bool isWhiteSpace() const;
    bool isLineTerminator() const;
    bool isHexDigit(ushort c) const;
    bool isOctalDigit(ushort c) const;

    int matchPunctuator(ushort c1, ushort c2,
                         ushort c3, ushort c4);
    ushort singleEscape(ushort c) const;
    ushort convertOctal(ushort c1, ushort c2,
                         ushort c3) const;
public:
    static unsigned char convertHex(ushort c1);
    static unsigned char convertHex(ushort c1, ushort c2);
    static QChar convertUnicode(ushort c1, ushort c2,
                                 ushort c3, ushort c4);
    static bool isIdentLetter(ushort c);
    static bool isDecimalDigit(ushort c);

    inline int ival() const { return qsyylval.ival; }
    inline double dval() const { return qsyylval.dval; }
    inline QScriptNameIdImpl *ustr() const { return qsyylval.ustr; }

    const QChar *characterBuffer() const { return buffer16; }
    int characterCount() const { return pos16; }

private:
    void record8(ushort c);
    void record16(QChar c);
    void recordStartPos();

    int findReservedWord(const QChar *buffer, int size) const;

    void syncProhibitAutomaticSemicolon();

    const QChar *code;
    uint length;
    int yycolumn;
    int startlineno;
    int startcolumn;
    int bol;     // begin of line

    union {
        int ival;
        double dval;
        QScriptNameIdImpl *ustr;
    } qsyylval;

    // current and following unicode characters
    ushort current, next1, next2, next3;

    struct keyword {
        const char *name;
        int token;
    };

    QString errmsg;
    Error err;

    bool wantRx;
    bool check_reserved;

    ParenthesesState parenthesesState;
    int parenthesesCount;
    bool prohibitAutomaticSemicolon;
};

} // namespace QScript

QT_END_NAMESPACE

#endif
