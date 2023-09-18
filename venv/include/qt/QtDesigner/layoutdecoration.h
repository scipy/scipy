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

#ifndef LAYOUTDECORATION_H
#define LAYOUTDECORATION_H

#include <QtDesigner/extension.h>

#include <QtCore/qobject.h>
#include <QtCore/qpair.h>

QT_BEGIN_NAMESPACE

class QPoint;
class QLayoutItem;
class QWidget;
class QRect;
class QLayout;

class QDesignerLayoutDecorationExtension
{
public:
    enum InsertMode
    {
        InsertWidgetMode,
        InsertRowMode,
        InsertColumnMode
    };

    virtual ~QDesignerLayoutDecorationExtension() {}

    virtual QList<QWidget*> widgets(QLayout *layout) const = 0;

    virtual QRect itemInfo(int index) const = 0;
    virtual int indexOf(QWidget *widget) const = 0;
    virtual int indexOf(QLayoutItem *item) const = 0;

    virtual InsertMode currentInsertMode() const = 0;
    virtual int currentIndex() const = 0;
    virtual QPair<int, int> currentCell() const = 0;
    virtual void insertWidget(QWidget *widget, const QPair<int, int> &cell) = 0;
    virtual void removeWidget(QWidget *widget) = 0;

    virtual void insertRow(int row) = 0;
    virtual void insertColumn(int column) = 0;
    virtual void simplify() = 0;

    virtual int findItemAt(const QPoint &pos) const = 0;
    virtual int findItemAt(int row, int column) const = 0; // atm only for grid.

    virtual void adjustIndicator(const QPoint &pos, int index) = 0;
};
Q_DECLARE_EXTENSION_INTERFACE(QDesignerLayoutDecorationExtension, "org.qt-project.Qt.Designer.LayoutDecoration")

QT_END_NAMESPACE

#endif // LAYOUTDECORATION_H
