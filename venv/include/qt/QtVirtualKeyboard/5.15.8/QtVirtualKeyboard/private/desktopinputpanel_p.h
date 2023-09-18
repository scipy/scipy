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

#ifndef DESKTOPINPUTPANEL_P_H
#define DESKTOPINPUTPANEL_P_H

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

#include <QtVirtualKeyboard/private/appinputpanel_p.h>

QT_BEGIN_NAMESPACE

class QWindow;

namespace QtVirtualKeyboard {

class DesktopInputPanelPrivate;

class QVIRTUALKEYBOARD_EXPORT DesktopInputPanel : public AppInputPanel
{
    Q_OBJECT
    Q_DECLARE_PRIVATE(DesktopInputPanel)
public:
    explicit DesktopInputPanel(QObject *parent = nullptr);
    ~DesktopInputPanel();

    void show() override;
    void hide() override;
    bool isVisible() const override;

    void setInputRect(const QRect &inputRect) override;

public slots:
    void createView() override;
    void destroyView() override;

protected slots:
    void repositionView(const QRect &rect);
    void focusWindowChanged(QWindow *focusWindow);
    void focusWindowVisibleChanged(bool visible);
    void previewRectangleChanged();
    void previewVisibleChanged();

protected:
    void updateInputRegion();
};

} // namespace QtVirtualKeyboard
QT_END_NAMESPACE

#endif // DESKTOPINPUTPANEL_P_H
