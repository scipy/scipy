/****************************************************************************
**
** Copyright (C) 2018 The Qt Company Ltd.
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

#ifndef QVIRTUALKEYBOARDINPUTCONTEXT_P_H
#define QVIRTUALKEYBOARDINPUTCONTEXT_P_H

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

#include <QObject>
#include <QRectF>
#include <QInputMethodEvent>
#include <QtVirtualKeyboard/qvirtualkeyboardinputcontext.h>
#include <QtVirtualKeyboard/private/shadowinputcontext_p.h>

QT_BEGIN_NAMESPACE

namespace QtVirtualKeyboard {

// Boolean variable QT_VIRTUALKEYBOARD_FORCE_EVENTS_WITHOUT_FOCUS will enable virtual keyboard
// to send key events without having any text input in focus when an environment
// variable QT_VIRTUALKEYBOARD_FORCE_EVENTS_WITHOUT_FOCUS is set. This will also enable virtual
// keyboard shift key usage in the same no focus situation.
extern const bool QT_VIRTUALKEYBOARD_FORCE_EVENTS_WITHOUT_FOCUS;

class PlatformInputContext;
class ShiftHandler;
}
class QVirtualKeyboardInputEngine;
class QVirtualKeyboardInputContextPrivate;

class QVIRTUALKEYBOARD_EXPORT QVirtualKeyboardInputContextPrivate : public QObject
{
    Q_OBJECT
    Q_DECLARE_PUBLIC(QVirtualKeyboardInputContext)
    Q_PROPERTY(bool focus READ focus WRITE setFocus NOTIFY focusChanged)
    Q_PROPERTY(QRectF keyboardRectangle READ keyboardRectangle WRITE setKeyboardRectangle NOTIFY keyboardRectangleChanged)
    Q_PROPERTY(QRectF previewRectangle READ previewRectangle WRITE setPreviewRectangle NOTIFY previewRectangleChanged)
    Q_PROPERTY(bool previewVisible READ previewVisible WRITE setPreviewVisible NOTIFY previewVisibleChanged)
    Q_PROPERTY(QString locale READ locale WRITE setLocale NOTIFY localeChanged)
    Q_PROPERTY(QObject *inputItem READ inputItem NOTIFY inputItemChanged)
    Q_PROPERTY(QtVirtualKeyboard::ShiftHandler *shiftHandler READ shiftHandler CONSTANT)
    Q_PROPERTY(QtVirtualKeyboard::ShadowInputContext *shadow READ shadow CONSTANT)
    Q_PROPERTY(QStringList inputMethods READ inputMethods CONSTANT)

    explicit QVirtualKeyboardInputContextPrivate(QVirtualKeyboardInputContext *q_ptr);
    void init();

public:
    enum class State {
        Reselect = 0x1,
        InputMethodEvent = 0x2,
        KeyEvent = 0x4,
        InputMethodClick = 0x8,
        SyncShadowInput = 0x10
    };
    Q_FLAG(State)
    Q_DECLARE_FLAGS(StateFlags, QVirtualKeyboardInputContextPrivate::State)

    ~QVirtualKeyboardInputContextPrivate();

    bool focus() const;
    void setFocus(bool focus);
    QRectF keyboardRectangle() const;
    void setKeyboardRectangle(QRectF rectangle);
    QRectF previewRectangle() const;
    void setPreviewRectangle(QRectF rectangle);
    bool previewVisible() const;
    void setPreviewVisible(bool visible);
    QString locale() const;
    void setLocale(const QString &locale);
    QObject *inputItem() const;
    QtVirtualKeyboard::ShiftHandler *shiftHandler() const;
    QtVirtualKeyboard::ShadowInputContext *shadow() const;
    QStringList inputMethods() const;

    // Helper functions
    Q_INVOKABLE bool fileExists(const QUrl &fileUrl);
    Q_INVOKABLE bool hasEnterKeyAction(QObject *item) const;
    Q_INVOKABLE void registerInputPanel(QObject *inputPanel);

Q_SIGNALS:
    void focusChanged();
    void keyboardRectangleChanged();
    void previewRectangleChanged();
    void previewVisibleChanged();
    void localeChanged();
    void inputItemChanged();
    void navigationKeyPressed(int key, bool isAutoRepeat);
    void navigationKeyReleased(int key, bool isAutoRepeat);

public Q_SLOTS:
    void hideInputPanel();
    void updateAvailableLocales(const QStringList &availableLocales);
    void forceCursorPosition(int anchorPosition, int cursorPosition);

private Q_SLOTS:
    void onInputItemChanged();

private:
    void sendPreedit(const QString &text, const QList<QInputMethodEvent::Attribute> &attributes, int replaceFrom, int replaceLength);
    void sendInputMethodEvent(QInputMethodEvent *event);
    void reset();
    void commit();
    void update(Qt::InputMethodQueries queries);
    void invokeAction(QInputMethod::Action action, int cursorPosition);
    bool filterEvent(const QEvent *event);
    void addSelectionAttribute(QList<QInputMethodEvent::Attribute> &attributes);
    bool testAttribute(const QList<QInputMethodEvent::Attribute> &attributes, QInputMethodEvent::AttributeType attributeType) const;
    int findAttribute(const QList<QInputMethodEvent::Attribute> &attributes, QInputMethodEvent::AttributeType attributeType) const;
    inline void setState(const State &state) { stateFlags.setFlag(state); }
    inline void clearState(const State &state) { stateFlags &= ~StateFlags(state); }
    inline bool testState(const State &state) const { return stateFlags.testFlag(state); }
    inline bool isEmptyState() const { return !stateFlags; }

private:
    QVirtualKeyboardInputContext *q_ptr;
    QtVirtualKeyboard::PlatformInputContext *platformInputContext;
    QVirtualKeyboardInputEngine *inputEngine;
    QtVirtualKeyboard::ShiftHandler *_shiftHandler;
    QPointer<QObject> inputPanel;
    QPointer<QObject> inputPanelParentItem;
    qreal inputPanelZ = .0;
    QRectF keyboardRect;
    QRectF previewRect;
    bool _previewVisible;
    bool animating;
    bool _focus;
    StateFlags stateFlags;
    int cursorPosition;
    int anchorPosition;
    int forceAnchorPosition;
    int _forceCursorPosition;
    Qt::InputMethodHints inputMethodHints;
    QString preeditText;
    QList<QInputMethodEvent::Attribute> preeditTextAttributes;
    QString surroundingText;
    QString selectedText;
    QRectF anchorRectangle;
    QRectF cursorRectangle;
    bool selectionControlVisible;
    bool anchorRectIntersectsClipRect;
    bool cursorRectIntersectsClipRect;
#ifdef QT_VIRTUALKEYBOARD_ARROW_KEY_NAVIGATION
    QSet<int> activeNavigationKeys;
#endif
    QSet<quint32> activeKeys;
    QtVirtualKeyboard::ShadowInputContext _shadow;

    friend class QtVirtualKeyboard::PlatformInputContext;
    friend class QVirtualKeyboardScopedState;
};

class QVirtualKeyboardScopedState
{
    Q_DISABLE_COPY(QVirtualKeyboardScopedState)
public:
    QVirtualKeyboardScopedState(QVirtualKeyboardInputContextPrivate *d, QVirtualKeyboardInputContextPrivate::State state) :
        d(d),
        state(state)
    {
        d->setState(state);
    }

    ~QVirtualKeyboardScopedState()
    {
        d->clearState(state);
    }

private:
    QVirtualKeyboardInputContextPrivate *d;
    QVirtualKeyboardInputContextPrivate::State state;
};

Q_DECLARE_OPERATORS_FOR_FLAGS(QVirtualKeyboardInputContextPrivate::StateFlags)

QT_END_NAMESPACE

Q_DECLARE_METATYPE(QVirtualKeyboardInputContextPrivate::State)

#endif // QVIRTUALKEYBOARDINPUTCONTEXT_P_H
