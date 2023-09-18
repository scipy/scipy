/****************************************************************************
**
** Copyright (C) 2018 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtSCriptTools module of the Qt Toolkit.
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

#ifndef QSCRIPTEDIT_P_H
#define QSCRIPTEDIT_P_H

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

#include <QtWidgets/qplaintextedit.h>

#include <QtCore/qhash.h>
#include <QtCore/qset.h>

QT_BEGIN_NAMESPACE

class QScriptEditExtraArea;

class Q_AUTOTEST_EXPORT QScriptEdit : public QPlainTextEdit
{
    Q_OBJECT
public:
    QScriptEdit(QWidget *parent = 0);
    ~QScriptEdit();

    int baseLineNumber() const;
    void setBaseLineNumber(int base);

    int executionLineNumber() const;
    void setExecutionLineNumber(int lineNumber, bool error);
    void setExecutableLineNumbers(const QSet<int> &lineNumbers);
    bool isExecutableLine(int lineNumber) const;

    int currentLineNumber() const;
    void gotoLine(int lineNumber);

    void setBreakpoint(int lineNumber);
    void setBreakpointEnabled(int lineNumber, bool enable);
    void deleteBreakpoint(int lineNumber);

    int extraAreaWidth() const;

signals:
    void breakpointToggleRequest(int lineNumber, bool on);
    void breakpointEnableRequest(int lineNumber, bool enable);

protected:
    void paintEvent(QPaintEvent *e);
    void resizeEvent(QResizeEvent *e);

    void extraAreaPaintEvent(QPaintEvent *e);
    void extraAreaMouseEvent(QMouseEvent *e);
    bool extraAreaEvent(QEvent *e);

private slots:
    void updateExtraAreaWidth();
    void updateExtraArea(const QRect &, int);
    void highlightCurrentLine();

private:
    QTextEdit::ExtraSelection currentLineSelection() const;
    QTextEdit::ExtraSelection currentExecutionLineSelection() const;
    void updateExtraSelections();

private:
    QScriptEditExtraArea *m_extraArea;
    int m_baseLineNumber;
    int m_executionLineNumber;
    QSet<int> m_executableLineNumbers;
    bool m_executionLineNumberHasError;
    int m_extraAreaToggleBlockNumber;

    struct BreakpointData {
        BreakpointData() : enabled(true) {}
        bool enabled;
    };
    QHash<int, BreakpointData> m_breakpoints;

private:
    Q_DISABLE_COPY(QScriptEdit)
    friend class QScriptEditExtraArea;
};

QT_END_NAMESPACE

#endif
