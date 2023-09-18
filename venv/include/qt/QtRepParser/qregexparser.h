/****************************************************************************
** Copyright (C) 2017 Ford Motor Company.
** All rights reserved.
**
** Copyright (C) 2017 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtRemoteObjects module of the Qt Toolkit.
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

#ifndef QREGEXPARSER_H
#define QREGEXPARSER_H

#include <QtCore/qshareddata.h>
#include <QtCore/qvarlengtharray.h>
#include <QtCore/qvariant.h>
#ifdef QT_BOOTSTRAPPED
#  include <QtCore/qregexp.h>
#  define REGEX QRegExp
#else
#  include <QtCore/qregularexpression.h>
#  define REGEX QRegularExpression
#endif
#include <QtCore/qmap.h>
#include <QtCore/qfile.h>
#include <QtCore/qtextstream.h>
#include <QtCore/qdebug.h>

struct MatchCandidate {
    MatchCandidate(const QString &n, const QString &t, int i) : name(n), matchText(t), index(i) {}
    QString name;
    QString matchText;
    int index;
};

QT_BEGIN_NAMESPACE

template <typename _Parser, typename _Table>
class QRegexParser: protected _Table
{
public:
    QRegexParser(int maxMatchLen=4096);
    virtual ~QRegexParser();

    virtual bool parse();

    virtual void reset() {}

    inline QVariant &sym(int index);

    void setBuffer(const QString &buffer);

    void setBufferFromDevice(QIODevice *device);

    void setDebug();

    QString errorString() const
    {
        return m_errorString;
    }

    void setErrorString(const QString &error)
    {
        m_errorString = error;
        qWarning() << m_errorString;
    }

    inline const QMap<QString, QString>& captured() const
    {
        return m_captured;
    }

    inline bool isDebug() const
    {
        return m_debug;
    }

    inline int lineNumber() const
    {
        return m_lineno;
    }

private:
    int nextToken();

    inline bool consumeRule(int rule)
    {
        return static_cast<_Parser*> (this)->consumeRule(rule);
    }

    enum { DefaultStackSize = 128 };

    struct Data: public QSharedData
    {
        Data(): stackSize (DefaultStackSize), tos (0) {}

        QVarLengthArray<int, DefaultStackSize> stateStack;
        QVarLengthArray<QVariant, DefaultStackSize> parseStack;
        int stackSize;
        int tos;

        void reallocateStack() {
            stackSize <<= 1;
            stateStack.resize(stackSize);
            parseStack.resize(stackSize);
        }
    };

    inline QString escapeString(QString s)
    {
        return s.replace(QLatin1Char('\n'), QLatin1String("\\n")).replace(QLatin1Char('\t'), QLatin1String("\\t"));
    }

    QSharedDataPointer<Data> d;

    QList<REGEX> m_regexes;
#ifndef QT_BOOTSTRAPPED
    QMap<QChar, QList<int> > regexCandidates;
#endif
    QList<int> m_tokens;
    QString m_buffer, m_lastMatchText;
    int m_loc, m_lastNewlinePosition;
    int m_lineno;
    int m_debug;
    QStringList m_tokenNames;
    QMap<QString, QString> m_captured;
    int m_maxMatchLen;
    QString m_errorString;
    QVector<QMap<int, QString> > m_names; //storage for match names
};

template <typename _Parser, typename _Table>
inline QVariant &QRegexParser<_Parser, _Table>::sym(int n)
{
    return d->parseStack [d->tos + n - 1];
}

template <typename _Parser, typename _Table>
QRegexParser<_Parser, _Table>::~QRegexParser()
{
}

template <typename _Parser, typename _Table>
bool QRegexParser<_Parser, _Table>::parse()
{
    m_errorString.clear();
    reset();
    const int INITIAL_STATE = 0;

    d->tos = 0;
    d->reallocateStack();

    int act = d->stateStack[++d->tos] = INITIAL_STATE;
    int token = -1;

    Q_FOREVER {
        if (token == -1 && - _Table::TERMINAL_COUNT != _Table::action_index[act])
            token = nextToken();

        act = _Table::t_action(act, token);

        if (d->stateStack[d->tos] == _Table::ACCEPT_STATE)
            return true;

        else if (act > 0) {
            if (++d->tos == d->stackSize)
                d->reallocateStack();

            d->parseStack[d->tos] = d->parseStack[d->tos - 1];
            d->stateStack[d->tos] = act;
            token = -1;
        }

        else if (act < 0) {
            int r = - act - 1;
            d->tos -= _Table::rhs[r];
            act = d->stateStack[d->tos++];
            if (!consumeRule(r))
                return false;
            act = d->stateStack[d->tos] = _Table::nt_action(act, _Table::lhs[r] - _Table::TERMINAL_COUNT);
        }

        else break;
    }

    setErrorString(QStringLiteral("Unknown token encountered"));
    return false;
}

template <typename _Parser, typename _Table>
QRegexParser<_Parser, _Table>::QRegexParser(int maxMatchLen) : d(new Data()), m_loc(0), m_lastNewlinePosition(0), m_lineno(1), m_debug(0), m_maxMatchLen(maxMatchLen)
{
    REGEX re(QStringLiteral("\\[([_a-zA-Z][_0-9a-zA-Z]*)(,\\s*M)?\\](.+)$"));
#ifdef QT_BOOTSTRAPPED
    REGEX nameMatch(QStringLiteral("\\((\\?<(.*)>).+\\)"));
    nameMatch.setMinimal(true);
#else
    re.optimize();
#endif
    QMap<QString, int> token_lookup;
    QMap<int, QString> names;
    for (int i = 1; i < _Table::lhs[0]; i++) {
        const QString text = QLatin1String(_Table::spell[i]);
        names.clear();
#ifdef QT_BOOTSTRAPPED
        if (re.indexIn(text) == 0) {
            const QString token = re.cap(1);
            const bool multiline = re.cap(2).length() > 0;
            QString pattern = re.cap(3);
            //We need to identify/remove any match names in the pattern, since
            //QRegExp doesn't support that feature
            int pos = 0, counter = 1, loc = nameMatch.indexIn(pattern, pos);
            while (loc >= 0) {
                const QString res = nameMatch.cap(2);
                if (!res.isEmpty()) {
                    names.insert(counter, res);
                    pattern.remove(nameMatch.cap(1));
                }
                pos += loc + nameMatch.matchedLength() - nameMatch.cap(1).length();
                loc = nameMatch.indexIn(pattern, pos);
                ++counter;
            }
            //We need to use indexIn, but that will search past the location we
            //pass in.  So prepend '^' and use QRegExp::CaretAtOffset.
            if (pattern.at(0) != QChar(QLatin1Char('^')))
                pattern.prepend(QChar(QLatin1Char('^')));
#else
        QRegularExpressionMatch match = re.match(text, 0, QRegularExpression::NormalMatch, QRegularExpression::DontCheckSubjectStringMatchOption);
        if (match.hasMatch()) {
            const QString token = match.captured(1);
            const bool multiline = match.captured(2).length() > 0;
            const QString pattern = match.captured(3);
#endif
            m_tokenNames.append(token);
            int index = i;
            if (token_lookup.contains(token))
                index = token_lookup[token];
            else
                token_lookup[token] = i;
#ifdef QT_BOOTSTRAPPED
            if (multiline)
                qWarning() << "The multiline grammar option is ignore in force_bootstrap mode.";
#endif
            REGEX pat(pattern);
#ifndef QT_BOOTSTRAPPED
            if (multiline)
                pat.setPatternOptions(QRegularExpression::DotMatchesEverythingOption);
#endif
            if (!pat.isValid())
                qCritical() << "Pattern error for token #" << i << "for" << text << "pattern =" << pat << ":" << pat.errorString();
            else {
#ifndef QT_BOOTSTRAPPED
                pat.optimize();
                int counter = 0;
                const auto namedCaptureGroups = pat.namedCaptureGroups();
                for (const QString &name : namedCaptureGroups) {
                    if (!name.isEmpty())
                        names.insert(counter, name);
                    ++counter;
                }
#endif
                m_names.append(names);
                m_regexes.append(pat);
                if (token.startsWith(QLatin1String("ignore")))
                    m_tokens.append(-1);
                else
                    m_tokens.append(index);
            }
        } else {
            qCritical() << "Error parsing regex at token #" << i << "for" << text << "Invalid syntax";
        }
    }
}

template <typename _Parser, typename _Table>
void QRegexParser<_Parser, _Table>::setBuffer(const QString &buffer)
{
    m_buffer = buffer;
}

template <typename _Parser, typename _Table>
void QRegexParser<_Parser, _Table>::setBufferFromDevice(QIODevice *device)
{
    QTextStream in(device);
    m_buffer = in.readAll();
}

template <typename _Parser, typename _Table>
void QRegexParser<_Parser, _Table>::setDebug()
{
    m_debug = true;
    for (int r = 0; r < _Table::RULE_COUNT; ++r)
    {
        int ridx = _Table::rule_index[r];
        int _rhs = _Table::rhs[r];
        qDebug("%3d) %s ::=", r + 1, _Table::spell[_Table::rule_info[ridx]]);
        ++ridx;
        for (int i = ridx; i < ridx + _rhs; ++i)
        {
            int symbol = _Table::rule_info[i];
            if (symbol > 0 && symbol < _Table::lhs[0])
                qDebug("     token_%s (pattern = %s)",qPrintable(m_tokenNames[symbol-1]),qPrintable(m_regexes[symbol-1].pattern()));
            else if (const char *name = _Table::spell[symbol])
                qDebug("     %s", name);
            else
                qDebug("     #%d", symbol);
        }
        qDebug();
    }
}

template <typename _Parser, typename _Table>
int QRegexParser<_Parser, _Table>::nextToken()
{
    static const REGEX newline(QLatin1String("(\\n)"));
    int token = -1;
    while (token < 0)
    {
        if (m_loc == m_buffer.size())
            return _Table::EOF_SYMBOL;

        //Check m_lastMatchText for newlines and update m_lineno
        //This isn't necessary, but being able to provide the line # and character #
        //where the match is failing sure makes building/debugging grammars easier.
#ifdef QT_BOOTSTRAPPED
        int loc = 0, pos = newline.indexIn(m_lastMatchText, loc);
        while (pos >= 0) {
            m_lineno++;
            loc += pos + 1;
            m_lastNewlinePosition += pos + 1;
            pos = newline.indexIn(m_lastMatchText, loc);
        }
#else //QT_BOOTSTRAPPED
        QRegularExpressionMatchIterator  matches = newline.globalMatch(m_lastMatchText);
        while (matches.hasNext()) {
            m_lineno++;
            QRegularExpressionMatch match = matches.next();
            if (!matches.hasNext())
                m_lastNewlinePosition += match.capturedEnd();
        }
#endif //!QT_BOOTSTRAPPED
        if (m_debug) {
            qDebug();
            qDebug() << "nextToken loop, line =" << m_lineno
                << "line position =" << m_loc - m_lastNewlinePosition
                << "next 5 characters =" << escapeString(m_buffer.mid(m_loc, 5));
        }
        int best = -1, maxLen = -1;
#ifndef QT_BOOTSTRAPPED
        QRegularExpressionMatch bestRegex;
#endif

        //Find the longest match.
        //If more than one are the same (longest) length, return the first one in
        //the order defined.
        QList<MatchCandidate> candidates;
#ifndef QT_BOOTSTRAPPED
        {
            //We used PCRE's PartialMatch to eliminate most of the regexes by the first
            //character, so we keep a regexCandidates map with the list of possible regexes
            //based on initial characters found so far.
            const QChar nextChar = m_buffer.at(m_loc);
            //Populate the list if we haven't seeen this character before
            if (!regexCandidates.contains(nextChar)) {
#  if (QT_VERSION >= QT_VERSION_CHECK(5, 5, 0))
                const QStringRef tmp = m_buffer.midRef(m_loc,1);
#  else
                const QString tmp = m_buffer.mid(m_loc,1);
#  endif
                int i = 0;
                regexCandidates[nextChar] = QList<int>();
                for (const QRegularExpression &re : qAsConst(m_regexes))
                {
                    QRegularExpressionMatch match = re.match(tmp, 0, QRegularExpression::PartialPreferFirstMatch, QRegularExpression::DontCheckSubjectStringMatchOption);
                    //qDebug() << nextChar << tmp << match.hasMatch() << match.hasPartialMatch() << re.pattern();
                    if (match.hasMatch() || match.hasPartialMatch())
                        regexCandidates[nextChar] << i;
                    i++;
                }
            }
            const auto indices = regexCandidates.value(nextChar);
            for (int i : indices)
            {
                //Seems like I should be able to run the regex on the entire string, but performance is horrible
                //unless I use a substring.
                //QRegularExpressionMatch match = m_regexes[i].match(m_buffer, m_loc, QRegularExpression::NormalMatch, QRegularExpression::AnchoredMatchOption);
#  if (QT_VERSION >= QT_VERSION_CHECK(5, 5, 0))
                QRegularExpressionMatch match = m_regexes.at(i).match(m_buffer.midRef(m_loc, m_maxMatchLen), 0, QRegularExpression::NormalMatch, QRegularExpression::AnchoredMatchOption | QRegularExpression::DontCheckSubjectStringMatchOption);
#  else
                QRegularExpressionMatch match = m_regexes.at(i).match(m_buffer.mid(m_loc, m_maxMatchLen), 0, QRegularExpression::NormalMatch, QRegularExpression::AnchoredMatchOption | QRegularExpression::DontCheckSubjectStringMatchOption);
#  endif
                if (match.hasMatch()) {
                    if (m_debug)
                        candidates << MatchCandidate(m_tokenNames[i], match.captured(), i);
                    if (match.capturedLength() > maxLen) {
                        best = i;
                        maxLen = match.capturedLength();
                        bestRegex = match;
                    }
                }
            }
        }
#else
        {
            int i = 0;
            for (const QRegExp &r : qAsConst(m_regexes))
            {
                if (r.indexIn(m_buffer, m_loc, QRegExp::CaretAtOffset) == m_loc) {
                    if (m_debug)
                        candidates << MatchCandidate(m_tokenNames[i], r.cap(0), i);
                    if (r.matchedLength() > maxLen) {
                        best = i;
                        maxLen = r.matchedLength();
                    }
                }
                ++i;
            }
        }
#endif
        if (best < 0) {
            setErrorString(QLatin1String("Error generating tokens from file, next characters >%1<").arg(m_buffer.midRef(m_loc, 15)));
            return -1;
        } else {
            const QMap<int, QString> &map = m_names.at(best);
            if (!map.isEmpty())
                m_captured.clear();
            for (auto iter = map.cbegin(), end = map.cend(); iter != end; ++iter) {
#ifdef QT_BOOTSTRAPPED
                m_captured.insert(iter.value(), m_regexes.at(best).cap(iter.key()));
#else
                m_captured.insert(iter.value(), bestRegex.captured(iter.key()));
#endif
            }
            if (m_debug) {
                qDebug() << "Match candidates:";
                for (const MatchCandidate &m : qAsConst(candidates)) {
                    QLatin1String result = m.index == best ? QLatin1String(" * ") : QLatin1String("   ");
                    qDebug() << qPrintable(result) << qPrintable(m.name) << qPrintable(escapeString(m.matchText));
                }
            }
            m_loc += maxLen;
            if (m_tokens.at(best) >= 0)
                token = m_tokens.at(best);
#ifdef QT_BOOTSTRAPPED
            m_lastMatchText = m_regexes.at(best).cap(0);
#else
            m_lastMatchText = bestRegex.captured(0);
#endif
        }
    }
    return token;
}

QT_END_NAMESPACE

#endif // QREGEXPARSER_H
