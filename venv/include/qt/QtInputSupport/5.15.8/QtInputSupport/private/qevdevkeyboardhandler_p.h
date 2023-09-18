/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtGui module of the Qt Toolkit.
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

#ifndef QEVDEVKEYBOARDHANDLER_P_H
#define QEVDEVKEYBOARDHANDLER_P_H

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

#include <qobject.h>
#include <QTimer>
#include <QDataStream>

#include <memory>

QT_BEGIN_NAMESPACE

class QSocketNotifier;

namespace QEvdevKeyboardMap {
    const quint32 FileMagic = 0x514d4150; // 'QMAP'

    struct Mapping {
        quint16 keycode;
        quint16 unicode;
        quint32 qtcode;
        quint8 modifiers;
        quint8 flags;
        quint16 special;

    };

    enum Flags {
        IsDead     = 0x01,
        IsLetter   = 0x02,
        IsModifier = 0x04,
        IsSystem   = 0x08
    };

    enum System {
        SystemConsoleFirst    = 0x0100,
        SystemConsoleMask     = 0x007f,
        SystemConsoleLast     = 0x017f,
        SystemConsolePrevious = 0x0180,
        SystemConsoleNext     = 0x0181,
        SystemReboot          = 0x0200,
        SystemZap             = 0x0300
    };

    struct Composing {
        quint16 first;
        quint16 second;
        quint16 result;
    };

    enum Modifiers {
        ModPlain   = 0x00,
        ModShift   = 0x01,
        ModAltGr   = 0x02,
        ModControl = 0x04,
        ModAlt     = 0x08,
        ModShiftL  = 0x10,
        ModShiftR  = 0x20,
        ModCtrlL   = 0x40,
        ModCtrlR   = 0x80
        // ModCapsShift = 0x100, // not supported!
    };
}

inline QDataStream &operator>>(QDataStream &ds, QEvdevKeyboardMap::Mapping &m)
{
    return ds >> m.keycode >> m.unicode >> m.qtcode >> m.modifiers >> m.flags >> m.special;
}

inline QDataStream &operator<<(QDataStream &ds, const QEvdevKeyboardMap::Mapping &m)
{
    return ds << m.keycode << m.unicode << m.qtcode << m.modifiers << m.flags << m.special;
}

inline QDataStream &operator>>(QDataStream &ds, QEvdevKeyboardMap::Composing &c)
{
    return ds >> c.first >> c.second >> c.result;
}

inline QDataStream &operator<<(QDataStream &ds, const QEvdevKeyboardMap::Composing &c)
{
    return ds << c.first << c.second << c.result;
}

class QFdContainer
{
    int m_fd;
    Q_DISABLE_COPY_MOVE(QFdContainer);
public:
    explicit QFdContainer(int fd = -1) noexcept : m_fd(fd) {}
    ~QFdContainer() { reset(); }

    int get() const noexcept { return m_fd; }

    int release() noexcept { int result = m_fd; m_fd = -1; return result; }
    void reset() noexcept;
    void reset(int fd) { reset(); m_fd = fd; }
};

class QEvdevKeyboardHandler : public QObject
{
public:
    QEvdevKeyboardHandler(const QString &device, QFdContainer &fd, bool disableZap, bool enableCompose, const QString &keymapFile);
    ~QEvdevKeyboardHandler();

    enum KeycodeAction {
        None               = 0,

        CapsLockOff        = 0x01000000,
        CapsLockOn         = 0x01000001,
        NumLockOff         = 0x02000000,
        NumLockOn          = 0x02000001,
        ScrollLockOff      = 0x03000000,
        ScrollLockOn       = 0x03000001,

        Reboot             = 0x04000000,

        PreviousConsole    = 0x05000000,
        NextConsole        = 0x05000001,
        SwitchConsoleFirst = 0x06000000,
        SwitchConsoleLast  = 0x0600007f,
        SwitchConsoleMask  = 0x0000007f
    };

    static std::unique_ptr<QEvdevKeyboardHandler> create(const QString &device,
                                         const QString &specification,
                                         const QString &defaultKeymapFile = QString());

    static Qt::KeyboardModifiers toQtModifiers(quint8 mod)
    {
        Qt::KeyboardModifiers qtmod = Qt::NoModifier;

        if (mod & (QEvdevKeyboardMap::ModShift | QEvdevKeyboardMap::ModShiftL | QEvdevKeyboardMap::ModShiftR))
            qtmod |= Qt::ShiftModifier;
        if (mod & (QEvdevKeyboardMap::ModControl | QEvdevKeyboardMap::ModCtrlL | QEvdevKeyboardMap::ModCtrlR))
            qtmod |= Qt::ControlModifier;
        if (mod & QEvdevKeyboardMap::ModAlt)
            qtmod |= Qt::AltModifier;

        return qtmod;
    }

    bool loadKeymap(const QString &file);
    void unloadKeymap();

    void readKeycode();
    KeycodeAction processKeycode(quint16 keycode, bool pressed, bool autorepeat);

    void switchLang();

private:
    void processKeyEvent(int nativecode, int unicode, int qtcode,
                         Qt::KeyboardModifiers modifiers, bool isPress, bool autoRepeat);
    void switchLed(int, bool);

    QString m_device;
    QFdContainer m_fd;
    QSocketNotifier *m_notify;

    // keymap handling
    quint8 m_modifiers;
    quint8 m_locks[3];
    int m_composing;
    quint16 m_dead_unicode;
    quint8 m_langLock;

    bool m_no_zap;
    bool m_do_compose;

    const QEvdevKeyboardMap::Mapping *m_keymap;
    int m_keymap_size;
    const QEvdevKeyboardMap::Composing *m_keycompose;
    int m_keycompose_size;

    static const QEvdevKeyboardMap::Mapping s_keymap_default[];
    static const QEvdevKeyboardMap::Composing s_keycompose_default[];
};


QT_END_NAMESPACE

#endif // QEVDEVKEYBOARDHANDLER_P_H
