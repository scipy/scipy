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

#ifndef DESIGNERWIDGETITEM_H
#define DESIGNERWIDGETITEM_H

#include "shared_global_p.h"

#include <QtWidgets/qlayoutitem.h>
#include <QtCore/qobject.h>

QT_BEGIN_NAMESPACE

class QDesignerFormEditorInterface;

namespace qdesigner_internal {

// QDesignerWidgetItem: A Layout Item that is used for non-containerextension-
// containers (QFrame, etc) on Designer forms. It prevents its widget
// from being slammed to size 0 if the widget has no layout:
// Pre 4.5, this item ensured only that QWidgets and QFrames were not squeezed
// to size 0 since they have an invalid size hint when non-laid out.
// Since 4.5, the item is used for every  non-containerextension-container.
// In case the container has itself a layout, it merely tracks the minimum
// size. If the container has no layout and is not subject to some stretch
// factor, it will return the last valid size. The effect is that after
// breaking a layout on a container within a layout, it just maintains its
// last size and is not slammed to 0,0. In addition, it can be resized.
// The class keeps track of the containing layout by tracking widget reparent
// and destroyed slots as Designer will for example re-create grid layouts to
// shrink them.

class QDESIGNER_SHARED_EXPORT QDesignerWidgetItem : public QObject, public QWidgetItemV2  {
    Q_DISABLE_COPY_MOVE(QDesignerWidgetItem)
    Q_OBJECT
public:
    explicit QDesignerWidgetItem(const QLayout *containingLayout, QWidget *w, Qt::Orientations o = Qt::Horizontal|Qt::Vertical);

    const QLayout *containingLayout() const;

    inline QWidget *constWidget() const { return const_cast<QDesignerWidgetItem*>(this)->widget(); }

    QSize minimumSize() const override;
    QSize sizeHint()    const override;

    // Resize: Takes effect if the contained widget does not have a layout
    QSize nonLaidOutMinSize() const;
    void setNonLaidOutMinSize(const QSize &s);

    QSize nonLaidOutSizeHint() const;
    void setNonLaidOutSizeHint(const QSize &s);

    // Check whether a QDesignerWidgetItem should be installed
    static bool check(const QLayout *layout, QWidget *w, Qt::Orientations *ptrToOrientations = nullptr);

    // Register itself using QLayoutPrivate's widget item factory method hook
    static void install();
    static void deinstall();

    // Check for a non-container extension container
    static bool isContainer(const QDesignerFormEditorInterface *core, QWidget *w);

    static bool subjectToStretch(const QLayout *layout, QWidget *w);

    bool eventFilter(QObject * watched, QEvent * event) override;

private slots:
    void layoutChanged();

private:
    void expand(QSize *s) const;
    bool subjectToStretch() const;

    const Qt::Orientations m_orientations;
    mutable QSize m_nonLaidOutMinSize;
    mutable QSize m_nonLaidOutSizeHint;
    mutable const QLayout *m_cachedContainingLayout;
};

// Helper class that ensures QDesignerWidgetItem is installed while an
// instance is in scope.

class QDESIGNER_SHARED_EXPORT QDesignerWidgetItemInstaller {
    Q_DISABLE_COPY_MOVE(QDesignerWidgetItemInstaller)

public:
    QDesignerWidgetItemInstaller();
    ~QDesignerWidgetItemInstaller();

private:
    static int m_instanceCount;
};

} // namespace qdesigner_internal

QT_END_NAMESPACE

#endif
