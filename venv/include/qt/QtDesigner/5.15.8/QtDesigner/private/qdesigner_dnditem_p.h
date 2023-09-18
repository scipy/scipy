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

#ifndef QDESIGNER_DNDITEM_H
#define QDESIGNER_DNDITEM_H

#include "shared_global_p.h"
#include <QtDesigner/abstractdnditem.h>

#include <QtCore/qpoint.h>
#include <QtCore/qlist.h>
#include <QtCore/qmimedata.h>

QT_BEGIN_NAMESPACE

class QDrag;
class QImage;
class QDropEvent;

namespace qdesigner_internal {

class QDESIGNER_SHARED_EXPORT QDesignerDnDItem: public QDesignerDnDItemInterface
{
public:
    explicit QDesignerDnDItem(DropType type, QWidget *source = nullptr);
    ~QDesignerDnDItem() override;

    DomUI *domUi() const override;
    QWidget *decoration() const override;
    QWidget *widget() const override;
    QPoint hotSpot() const override;
    QWidget *source() const override;

    DropType type() const override;

protected:
    void setDomUi(DomUI *dom_ui);
    void init(DomUI *ui, QWidget *widget, QWidget *decoration, const QPoint &global_mouse_pos);

private:
    QWidget *m_source;
    const DropType m_type;
    const QPoint m_globalStartPos;
    DomUI *m_dom_ui;
    QWidget *m_widget;
    QWidget *m_decoration;
    QPoint m_hot_spot;

    Q_DISABLE_COPY_MOVE(QDesignerDnDItem)
};

// Mime data for use with designer drag and drop operations.

class  QDESIGNER_SHARED_EXPORT QDesignerMimeData : public QMimeData {
    Q_OBJECT

public:
    using QDesignerDnDItems = QList<QDesignerDnDItemInterface *>;

    ~QDesignerMimeData() override;

    const QDesignerDnDItems &items() const { return m_items; }

    // Execute a drag and drop operation.
    static Qt::DropAction execDrag(const QDesignerDnDItems &items, QWidget * dragSource);

    QPoint hotSpot() const { return m_hotSpot; }

    // Move the decoration. Required for drops over form windows as the position
    // is derived from the decoration position.
    void moveDecoration(const QPoint &globalPos) const;

    // For a move operation, create the undo command sequence to remove
    // the widgets from the source form.
    static void removeMovedWidgetsFromSourceForm(const QDesignerDnDItems &items);

    // Accept an event with the proper action.
    void acceptEvent(QDropEvent *e) const;

    // Helper to accept an event with the desired action.
    static void acceptEventWithAction(Qt::DropAction desiredAction, QDropEvent *e);

private:
    QDesignerMimeData(const QDesignerDnDItems &items, QDrag *drag);
    Qt::DropAction proposedDropAction() const;

    static void setImageTransparency(QImage &image, int alpha);

    const QDesignerDnDItems m_items;
    QPoint m_globalStartPos;
    QPoint m_hotSpot;
};

} // namespace qdesigner_internal

QT_END_NAMESPACE

#endif // QDESIGNER_DNDITEM_H
