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

#ifndef FORMLAYOUTMENU
#define FORMLAYOUTMENU

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

#include "shared_global_p.h"
#include <QtCore/qobject.h>
#include <QtCore/qlist.h>
#include <QtCore/qpointer.h>

QT_BEGIN_NAMESPACE

class QDesignerFormWindowInterface;

class QAction;
class QWidget;

namespace qdesigner_internal {

// Task menu to be used for form layouts. Offers an options "Add row" which
// pops up a dialog in which the user can specify label name, text and buddy.
class QDESIGNER_SHARED_EXPORT FormLayoutMenu : public QObject
{
    Q_DISABLE_COPY_MOVE(FormLayoutMenu)
    Q_OBJECT
public:
    using ActionList = QList<QAction *>;

    explicit FormLayoutMenu(QObject *parent);

    // Populate a list of actions with the form layout actions.
    void populate(QWidget *w, QDesignerFormWindowInterface *fw, ActionList &actions);
    // For implementing QDesignerTaskMenuExtension::preferredEditAction():
    // Return appropriate action for double clicking.
    QAction *preferredEditAction(QWidget *w, QDesignerFormWindowInterface *fw);

private slots:
    void slotAddRow();

private:
    QAction *m_separator1;
    QAction *m_populateFormAction;
    QAction *m_separator2;
    QPointer<QWidget> m_widget;
};
}  // namespace qdesigner_internal

QT_END_NAMESPACE

#endif // FORMLAYOUTMENU
