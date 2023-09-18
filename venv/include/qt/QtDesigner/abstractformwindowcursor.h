/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt Designer of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:GPL-EXCEPT$
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
** General Public License version 3 as published by the Free Software
** Foundation with exceptions as appearing in the file LICENSE.GPL3-EXCEPT
** included in the packaging of this file. Please review the following
** information to ensure the GNU General Public License requirements will
** be met: https://www.gnu.org/licenses/gpl-3.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef ABSTRACTFORMWINDOWCURSOR_H
#define ABSTRACTFORMWINDOWCURSOR_H

#include <QtDesigner/sdk_global.h>

QT_BEGIN_NAMESPACE

class QDesignerFormWindowInterface;
class QWidget;
class QVariant;
class QString;

class QDESIGNER_SDK_EXPORT QDesignerFormWindowCursorInterface
{
public:
    enum MoveOperation
    {
        NoMove,

        Start,
        End,
        Next,
        Prev,
        Left,
        Right,
        Up,
        Down
    };

    enum MoveMode
    {
        MoveAnchor,
        KeepAnchor
    };

public:
    virtual ~QDesignerFormWindowCursorInterface() {}

    virtual QDesignerFormWindowInterface *formWindow() const = 0;

    virtual bool movePosition(MoveOperation op, MoveMode mode = MoveAnchor) = 0;

    virtual int position() const = 0;
    virtual void setPosition(int pos, MoveMode mode = MoveAnchor) = 0;

    virtual QWidget *current() const = 0;

    virtual int widgetCount() const = 0;
    virtual QWidget *widget(int index) const = 0;

    virtual bool hasSelection() const = 0;
    virtual int selectedWidgetCount() const = 0;
    virtual QWidget *selectedWidget(int index) const = 0;

    virtual void setProperty(const QString &name, const QVariant &value) = 0;
    virtual void setWidgetProperty(QWidget *widget, const QString &name, const QVariant &value) = 0;
    virtual void resetWidgetProperty(QWidget *widget, const QString &name) = 0;

    bool isWidgetSelected(QWidget *widget) const;
};

QT_END_NAMESPACE

#endif // ABSTRACTFORMWINDOWCURSOR_H
