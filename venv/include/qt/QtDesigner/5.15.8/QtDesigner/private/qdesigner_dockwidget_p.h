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

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API.  It exists for the convenience
// of Qt Designer.  This header
// file may change from version to version without notice, or even be removed.
//
// We mean it.
//

#ifndef QDESIGNER_DOCKWIDGET_H
#define QDESIGNER_DOCKWIDGET_H

#include "shared_global_p.h"

#include <qdesigner_propertysheet_p.h>

#include <QtWidgets/qdockwidget.h>

QT_BEGIN_NAMESPACE

class QDesignerFormWindowInterface;

class QDESIGNER_SHARED_EXPORT QDesignerDockWidget: public QDockWidget
{
    Q_OBJECT
    Q_PROPERTY(Qt::DockWidgetArea dockWidgetArea READ dockWidgetArea WRITE setDockWidgetArea DESIGNABLE true STORED false)
    Q_PROPERTY(bool docked READ docked WRITE setDocked DESIGNABLE true STORED false)
public:
    QDesignerDockWidget(QWidget *parent = nullptr);
    ~QDesignerDockWidget() override;

    bool docked() const;
    void setDocked(bool b);

    Qt::DockWidgetArea dockWidgetArea() const;
    void setDockWidgetArea(Qt::DockWidgetArea dockWidgetArea);

    bool inMainWindow() const;

private:
    QDesignerFormWindowInterface *formWindow() const;
    QMainWindow *findMainWindow() const;
};

class QDESIGNER_SHARED_EXPORT QDockWidgetPropertySheet : public QDesignerPropertySheet
{
     Q_OBJECT
public:
    using QDesignerPropertySheet::QDesignerPropertySheet;

    bool isEnabled(int index) const override;
};

using QDockWidgetPropertySheetFactory =
    QDesignerPropertySheetFactory<QDockWidget, QDockWidgetPropertySheet>;

QT_END_NAMESPACE

#endif // QDESIGNER_DOCKWIDGET_H
