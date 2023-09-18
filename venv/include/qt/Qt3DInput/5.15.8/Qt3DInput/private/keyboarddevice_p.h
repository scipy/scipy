/****************************************************************************
**
** Copyright (C) 2014 Klaralvdalens Datakonsult AB (KDAB).
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt3D module of the Qt Toolkit.
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

#ifndef QT3DINPUT_INPUT_KEYBOARDDEVICE_P_H
#define QT3DINPUT_INPUT_KEYBOARDDEVICE_P_H

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API.  It exists for the convenience
// of other Qt classes.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.
//

#include <Qt3DInput/QKeyEvent>
#include <Qt3DCore/qnodeid.h>

#include <Qt3DInput/private/handle_types_p.h>
#include <Qt3DInput/private/qabstractphysicaldevicebackendnode_p.h>

QT_BEGIN_NAMESPACE

namespace Qt3DInput {

class QInputAspect;

namespace Input {

class InputHandler;

class Q_AUTOTEST_EXPORT KeyboardDevice : public Qt3DInput::QAbstractPhysicalDeviceBackendNode
{
public:
    KeyboardDevice();
    void cleanup() override;

    void requestFocusForInput(Qt3DCore::QNodeId inputId);
    void setInputHandler(InputHandler *handler);

    void setCurrentFocusItem(Qt3DCore::QNodeId input);

    float axisValue(int axisIdentifier) const override;
    bool isButtonPressed(int buttonIdentifier) const override;

    void updateKeyEvents(const QList<QT_PREPEND_NAMESPACE(QKeyEvent)> &events);

    inline Qt3DCore::QNodeId currentFocusItem() const { return m_currentFocusItem; }
    inline Qt3DCore::QNodeId lastKeyboardInputRequester() const { return m_lastRequester; }

private:
    void setButtonValue(int key, bool value);

    InputHandler *m_inputHandler;
    QVector<Qt3DCore::QNodeId> m_keyboardInputs;
    Qt3DCore::QNodeId m_lastRequester;
    Qt3DCore::QNodeId m_currentFocusItem;

    union KeyStates {

        struct Buttons
        {
            // first 4 bytes
            bool keyEscape:1;           // 0
            bool keyTab:1;              // 1
            bool keyBacktab:1;          // 2
            bool keyBackspace:1;        // 3
            bool keyReturn:1;           // 4
            bool keyEnter:1;            // 5
            bool keyInsert:1;           // 6
            bool keyDelete:1;           // 7
            bool keyPause:1;            // 8
            bool keyPrint:1;            // 9
            bool keySysReq:1;           // 10
            bool keyClear:1;            // 11
            bool keyHome:1;             // 12
            bool keyEnd:1;              // 13
            bool keyLeft:1;             // 14
            bool keyUp:1;               // 15
            bool keyRight:1;            // 16
            bool keyDown:1;             // 17
            bool keyPageUp:1;           // 18
            bool keyPageDown:1;         // 19
            bool keyShift:1;            // 20
            bool keyControl:1;          // 21
            bool keyMeta:1;             // 22
            bool keyAlt:1;              // 23
            bool keyCapsLock:1;         // 24
            bool keyNumLock:1;          // 25
            bool keyScrollLock:1;       // 26
            bool keyF1:1;               // 27
            bool keyF2:1;               // 28
            bool keyF3:1;               // 29
            bool keyF4:1;               // 30
            bool keyF5:1;               // 31

            // second 4 bytes
            bool keyF6:1;               // 0
            bool keyF7:1;               // 1
            bool keyF8:1;               // 2
            bool keyF9:1;               // 3
            bool keyF10:1;              // 4
            bool keyF11:1;              // 5
            bool keyF12:1;              // 6
            bool keyF13:1;              // 7
            bool keyF14:1;              // 8
            bool keyF15:1;              // 9
            bool keyF16:1;              // 10
            bool keyF17:1;              // 11
            bool keyF18:1;              // 12
            bool keyF19:1;              // 13
            bool keyF20:1;              // 14
            bool keyF21:1;              // 15
            bool keyF22:1;              // 16
            bool keyF23:1;              // 17
            bool keyF24:1;              // 18
            bool keyF25:1;              // 19
            bool keyF26:1;              // 20
            bool keyF27:1;              // 21
            bool keyF28:1;              // 22
            bool keyF29:1;              // 23
            bool keyF30:1;              // 24
            bool keyF31:1;              // 25
            bool keyF32:1;              // 26
            bool keyF33:1;              // 27
            bool keyF34:1;              // 28
            bool keyF35:1;              // 29
            bool keySuper_L:1;          // 30
            bool keySuper_R:1;          // 31

            // third 4 bytes
            // unused                   // 0
            bool keyMenu:1;             // 1
            bool keyHyper_L:1;          // 2
            bool keyHyper_R:1;          // 3
            bool keyHelp:1;             // 4
            bool keyDirection_L:1;      // 5
            bool keyDirection_R:1;      // 6
            bool keySpace:1;            // 7
            bool keyAny:1;              // 8
            bool keyExclam:1;           // 9
            bool keyQuoteDbl:1;         // 10
            bool keyNumberSign:1;       // 11
            bool keyDollar:1;           // 12
            bool keyPercent:1;          // 13
            bool keyAmpersand:1;        // 14
            bool keyApostrophe:1;       // 15
            bool keyParenLeft:1;        // 16
            bool keyParenRight:1;       // 17
            bool keyAsterisk:1;         // 18
            bool keyPlus:1;             // 19
            bool keyComma:1;            // 20
            bool keyMinus:1;            // 21
            bool keyPeriod:1;           // 22
            bool keySlash:1;            // 23
            bool key0:1;                // 24
            bool key1:1;                // 25
            bool key2:1;                // 26
            bool key3:1;                // 27
            bool key4:1;                // 28
            bool key5:1;                // 29
            bool key6:1;                // 30
            bool key7:1;                // 31

            // fourth 4 bytes
            bool key8:1;                // 0
            bool key9:1;                // 1
            bool keyColon:1;            // 2
            bool keySemicolon:1;        // 3
            bool keyLess:1;             // 4
            bool keyEqual:1;            // 5
            bool keyGreater:1;          // 6
            bool keyQuestion:1;         // 7
            bool keyAt:1;               // 8
            bool keyA:1;                // 9
            bool keyB:1;                // 10
            bool keyC:1;                // 11
            bool keyD:1;                // 12
            bool keyE:1;                // 13
            bool keyF:1;                // 14
            bool keyG:1;                // 15
            bool keyH:1;                // 16
            bool keyI:1;                // 17
            bool keyJ:1;                // 18
            bool keyK:1;                // 19
            bool keyL:1;                // 20
            bool keyM:1;                // 21
            bool keyN:1;                // 22
            bool keyO:1;                // 23
            bool keyP:1;                // 24
            bool keyQ:1;                // 25
            bool keyR:1;                // 26
            bool keyS:1;                // 27
            bool keyT:1;                // 28
            bool keyU:1;                // 29
            bool keyV:1;                // 30
            bool keyW:1;                // 31

            // fifth 4 bytes
            bool keyX:1;                // 0
            bool keyY:1;                // 1
            bool keyZ:1;                // 2
            bool keyBracketLeft:1;      // 3
            bool keyBackslash:1;        // 4
            bool keyBracketRight:1;     // 5
            bool keyAsciiCircum:1;      // 6
            bool keyUnderscore:1;       // 7
            bool keyQuoteLeft:1;        // 8
            bool keyBraceLeft:1;        // 9
            bool keyBar:1;              // 10
            bool keyBraceRight:1;       // 11
            bool keyAsciiTilde:1;       // 12
            bool keyplusminus:1;        // 13
            bool keyonesuperior:1;      // 14
            bool keymultiply:1;         // 15
            bool keydivision:1;         // 16
            bool keyydiaeresis:1;       // 17
        };
        qint32 keys[5];
    };

    KeyStates m_keyStates;
};

class KeyboardDeviceFunctor : public Qt3DCore::QBackendNodeMapper
{
public:
    explicit KeyboardDeviceFunctor(QInputAspect *inputaspect, InputHandler *handler);

    Qt3DCore::QBackendNode *create(const Qt3DCore::QNodeCreatedChangeBasePtr &change) const override;
    Qt3DCore::QBackendNode *get(Qt3DCore::QNodeId id) const override;
    void destroy(Qt3DCore::QNodeId id) const override;

private:
    QInputAspect *m_inputAspect;
    InputHandler *m_handler;
};

} // namespace Input
} // namespace Qt3DInput

QT_END_NAMESPACE

#endif // QT3DINPUT_INPUT_KEYBOARDDEVICE_P_H
