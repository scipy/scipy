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

#ifndef PLATFORMINPUTCONTEXT_P_H
#define PLATFORMINPUTCONTEXT_P_H

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

#include <qevent.h>
#include <qpa/qplatforminputcontext.h>
#include <QPointer>
#include <QLocale>
#include <QtVirtualKeyboard/private/inputselectionhandle_p.h>
#include <QtVirtualKeyboard/private/desktopinputselectioncontrol_p.h>

QT_BEGIN_NAMESPACE

class QVirtualKeyboardInputContext;
class QVirtualKeyboardInputContextPrivate;

namespace QtVirtualKeyboard {

class AbstractInputPanel;

class QVIRTUALKEYBOARD_EXPORT PlatformInputContext : public QPlatformInputContext
{
    Q_OBJECT
public:
    explicit PlatformInputContext();
    ~PlatformInputContext();

    bool isValid() const override;

    void reset() override;
    void commit() override;
    void update(Qt::InputMethodQueries queries) override;
    void invokeAction(QInputMethod::Action action, int cursorPosition) override;
    QRectF keyboardRect() const override;

    bool isAnimating() const override;

    void showInputPanel() override;
    void hideInputPanel() override;
    bool isInputPanelVisible() const override;

    QLocale locale() const override;
    void setLocale(QLocale locale);
    Qt::LayoutDirection inputDirection() const override;
    void setInputDirection(Qt::LayoutDirection direction);

    QObject *focusObject();
    void setFocusObject(QObject *object) override;

    QVirtualKeyboardInputContext *inputContext() const;

    bool eventFilter(QObject *object, QEvent *event) override;

    void setInputMethods(const QStringList &inputMethods);
    QStringList inputMethods() const;

signals:
    void focusObjectChanged();

protected:
    void sendEvent(QEvent *event);
    void sendKeyEvent(QKeyEvent *event);
    QVariant inputMethodQuery(Qt::InputMethodQuery query);
    void setInputContext(QVirtualKeyboardInputContext *context);

private slots:
    void keyboardRectangleChanged();
    void updateInputPanelVisible();

private:
    friend class ::QVirtualKeyboardInputContext;
    friend class ::QVirtualKeyboardInputContextPrivate;
    QPointer<QVirtualKeyboardInputContext> m_inputContext;
    QStringList m_inputMethods;
    QPointer<AbstractInputPanel> m_inputPanel;
    QPointer<DesktopInputSelectionControl> m_selectionControl;
    QPointer<QObject> m_focusObject;
    QLocale m_locale;
    Qt::LayoutDirection m_inputDirection;
    QEvent *m_filterEvent;
    bool m_visible;
    bool m_desktopModeDisabled;
};

} // namespace QtVirtualKeyboard
QT_END_NAMESPACE

#endif
