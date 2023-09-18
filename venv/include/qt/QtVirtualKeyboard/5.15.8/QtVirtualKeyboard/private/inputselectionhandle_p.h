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

#ifndef INPUTSELECTIONHANDLE_P_H
#define INPUTSELECTIONHANDLE_P_H

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

#include <QtGui/QRasterWindow>
#include <QtVirtualKeyboard/qvirtualkeyboard_global.h>

QT_BEGIN_NAMESPACE

class QWindow;

namespace QtVirtualKeyboard {

class DesktopInputSelectionControl;

class QVIRTUALKEYBOARD_EXPORT InputSelectionHandle : public QRasterWindow
{
    Q_OBJECT

public:
    InputSelectionHandle(DesktopInputSelectionControl *control, QWindow *eventWindow);
    void applyImage(const QSize &windowSize);

protected:
    void paintEvent(QPaintEvent *pe) override;
    bool event(QEvent *event) override;

private:
    DesktopInputSelectionControl *m_control;
    QWindow *m_eventWindow;
};

} // namespace QtVirtualKeyboard
QT_END_NAMESPACE

#endif // INPUTSELECTIONHANDLE_P_H
