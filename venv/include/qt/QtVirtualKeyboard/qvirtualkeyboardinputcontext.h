/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt Virtual Keyboard module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:GPL$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 3 or (at your option) any later version
** approved by the KDE Free Qt Foundation. The licenses are as published by
** the Free Software Foundation and appearing in the file LICENSE.GPL3
** included in the packaging of this file. Please review the following
** information to ensure the GNU General Public License requirements will
** be met: https://www.gnu.org/licenses/gpl-3.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef QVIRTUALKEYBOARDINPUTCONTEXT_H
#define QVIRTUALKEYBOARDINPUTCONTEXT_H

#include <QObject>
#include <QRectF>
#include <QLocale>
#include <QInputMethodEvent>
#include <QInputMethod>
#include <QtVirtualKeyboard/qvirtualkeyboard_global.h>

QT_BEGIN_NAMESPACE

namespace QtVirtualKeyboard {
class PlatformInputContext;
}
class QVirtualKeyboardInputEngine;
class QVirtualKeyboardInputContextPrivate;

class QVIRTUALKEYBOARD_EXPORT QVirtualKeyboardInputContext : public QObject
{
    Q_OBJECT
    Q_DISABLE_COPY(QVirtualKeyboardInputContext)
    Q_DECLARE_PRIVATE(QVirtualKeyboardInputContext)
    Q_PROPERTY(bool shift READ isShiftActive NOTIFY shiftActiveChanged)
    Q_PROPERTY(bool shiftActive READ isShiftActive NOTIFY shiftActiveChanged REVISION 4)
    Q_PROPERTY(bool capsLock READ isCapsLockActive NOTIFY capsLockActiveChanged)
    Q_PROPERTY(bool capsLockActive READ isCapsLockActive NOTIFY capsLockActiveChanged REVISION 4)
    Q_PROPERTY(bool uppercase READ isUppercase NOTIFY uppercaseChanged)
    Q_PROPERTY(int anchorPosition READ anchorPosition NOTIFY anchorPositionChanged)
    Q_PROPERTY(int cursorPosition READ cursorPosition NOTIFY cursorPositionChanged)
    Q_PROPERTY(Qt::InputMethodHints inputMethodHints READ inputMethodHints NOTIFY inputMethodHintsChanged)
    Q_PROPERTY(QString preeditText READ preeditText WRITE setPreeditText NOTIFY preeditTextChanged)
    Q_PROPERTY(QString surroundingText READ surroundingText NOTIFY surroundingTextChanged)
    Q_PROPERTY(QString selectedText READ selectedText NOTIFY selectedTextChanged)
    Q_PROPERTY(QRectF anchorRectangle READ anchorRectangle NOTIFY anchorRectangleChanged)
    Q_PROPERTY(QRectF cursorRectangle READ cursorRectangle NOTIFY cursorRectangleChanged)
    Q_PROPERTY(bool animating READ isAnimating WRITE setAnimating NOTIFY animatingChanged)
    Q_PROPERTY(QString locale READ locale NOTIFY localeChanged)
    Q_PROPERTY(QObject *inputItem READ inputItem NOTIFY inputItemChanged)
    Q_PROPERTY(QVirtualKeyboardInputEngine *inputEngine READ inputEngine CONSTANT)
    Q_PROPERTY(bool selectionControlVisible READ isSelectionControlVisible NOTIFY selectionControlVisibleChanged)
    Q_PROPERTY(bool anchorRectIntersectsClipRect READ anchorRectIntersectsClipRect NOTIFY anchorRectIntersectsClipRectChanged)
    Q_PROPERTY(bool cursorRectIntersectsClipRect READ cursorRectIntersectsClipRect NOTIFY cursorRectIntersectsClipRectChanged)
    Q_PROPERTY(QVirtualKeyboardInputContextPrivate *priv READ priv CONSTANT)

public:
    explicit QVirtualKeyboardInputContext(QObject *parent = nullptr);
    ~QVirtualKeyboardInputContext();

    bool isShiftActive() const;
    bool isCapsLockActive() const;
    bool isUppercase() const;
    int anchorPosition() const;
    int cursorPosition() const;
    Qt::InputMethodHints inputMethodHints() const;
    QString preeditText() const;
    void setPreeditText(const QString &text, QList<QInputMethodEvent::Attribute> attributes = QList<QInputMethodEvent::Attribute>(), int replaceFrom = 0, int replaceLength = 0);
    QList<QInputMethodEvent::Attribute> preeditTextAttributes() const;
    QString surroundingText() const;
    QString selectedText() const;
    QRectF anchorRectangle() const;
    QRectF cursorRectangle() const;
    bool isAnimating() const;
    void setAnimating(bool isAnimating);
    QString locale() const;
    QObject *inputItem() const;
    QVirtualKeyboardInputEngine *inputEngine() const;
    bool isSelectionControlVisible() const;
    bool anchorRectIntersectsClipRect() const;
    bool cursorRectIntersectsClipRect() const;
    QVirtualKeyboardInputContextPrivate *priv() const;

    Q_INVOKABLE void sendKeyClick(int key, const QString &text, int modifiers = 0);
    Q_INVOKABLE void commit();
    Q_INVOKABLE void commit(const QString &text, int replaceFrom = 0, int replaceLength = 0);
    Q_INVOKABLE void clear();

    // For selection handles
    Q_INVOKABLE void setSelectionOnFocusObject(const QPointF &anchorPos, const QPointF &cursorPos);

Q_SIGNALS:
    void preeditTextChanged();
    void inputMethodHintsChanged();
    void surroundingTextChanged();
    void selectedTextChanged();
    void anchorPositionChanged();
    void cursorPositionChanged();
    void anchorRectangleChanged();
    void cursorRectangleChanged();
    void shiftActiveChanged();
    void capsLockActiveChanged();
    void uppercaseChanged();
    void animatingChanged();
    void localeChanged();
    void inputItemChanged();
    void selectionControlVisibleChanged();
    void anchorRectIntersectsClipRectChanged();
    void cursorRectIntersectsClipRectChanged();

private:

    QScopedPointer<QVirtualKeyboardInputContextPrivate> d_ptr;
};

QT_END_NAMESPACE

#endif
