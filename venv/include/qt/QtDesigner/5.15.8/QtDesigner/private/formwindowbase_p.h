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

#ifndef FORMWINDOWBASE_H
#define FORMWINDOWBASE_H

#include "shared_global_p.h"

#include <QtDesigner/abstractformwindow.h>

#include <QtCore/qvariant.h>
#include <QtCore/qlist.h>

QT_BEGIN_NAMESPACE

class QDesignerDnDItemInterface;
class QMenu;
class QtResourceSet;
class QDesignerPropertySheet;

namespace qdesigner_internal {

class QEditorFormBuilder;
class DeviceProfile;
class Grid;

class DesignerPixmapCache;
class DesignerIconCache;
class FormWindowBasePrivate;

class QDESIGNER_SHARED_EXPORT FormWindowBase: public QDesignerFormWindowInterface
{
    Q_OBJECT
public:
    enum HighlightMode  { Restore, Highlight };

    explicit FormWindowBase(QDesignerFormEditorInterface *core, QWidget *parent = nullptr,
                            Qt::WindowFlags flags = {});
    ~FormWindowBase() override;

    QVariantMap formData();
    void setFormData(const QVariantMap &vm);

    QStringList checkContents() const override;

    // Deprecated
    QPoint grid() const override;

    // Deprecated
    void setGrid(const QPoint &grid) override;

    bool hasFeature(Feature f) const override;
    Feature features() const override;
    void setFeatures(Feature f) override;

    const Grid &designerGrid() const;
    void setDesignerGrid(const  Grid& grid);

    bool hasFormGrid() const;
    void setHasFormGrid(bool b);

    bool gridVisible() const;

    ResourceFileSaveMode resourceFileSaveMode() const override;
    void setResourceFileSaveMode(ResourceFileSaveMode behavior) override;

    static const Grid &defaultDesignerGrid();
    static void setDefaultDesignerGrid(const Grid& grid);

    // Overwrite to initialize and return a full popup menu for a managed widget
    virtual QMenu *initializePopupMenu(QWidget *managedWidget);
    // Helper to create a basic popup menu from task menu extensions (internal/public)
    static QMenu *createExtensionTaskMenu(QDesignerFormWindowInterface *fw, QObject *o, bool trailingSeparator = true);

    virtual bool dropWidgets(const QList<QDesignerDnDItemInterface*> &item_list, QWidget *target,
                             const QPoint &global_mouse_pos) = 0;

    // Helper to find the widget at the mouse position with some flags.
    enum WidgetUnderMouseMode { FindSingleSelectionDropTarget, FindMultiSelectionDropTarget };
    QWidget *widgetUnderMouse(const QPoint &formPos, WidgetUnderMouseMode m);

    virtual QWidget *widgetAt(const QPoint &pos) = 0;
    virtual QWidget *findContainer(QWidget *w, bool excludeLayout) const = 0;

    void deleteWidgetList(const QWidgetList &widget_list);

    virtual void highlightWidget(QWidget *w, const QPoint &pos, HighlightMode mode = Highlight) = 0;

    enum PasteMode { PasteAll, PasteActionsOnly };
#if QT_CONFIG(clipboard)
    virtual void paste(PasteMode pasteMode) = 0;
#endif

    // Factory method to create a form builder
    virtual QEditorFormBuilder *createFormBuilder() = 0;

    virtual bool blockSelectionChanged(bool blocked) = 0;

    DesignerPixmapCache *pixmapCache() const;
    DesignerIconCache *iconCache() const;
    QtResourceSet *resourceSet() const override;
    void setResourceSet(QtResourceSet *resourceSet) override;
    void addReloadableProperty(QDesignerPropertySheet *sheet, int index);
    void removeReloadableProperty(QDesignerPropertySheet *sheet, int index);
    void addReloadablePropertySheet(QDesignerPropertySheet *sheet, QObject *object);
    void reloadProperties();

    void emitWidgetRemoved(QWidget *w);
    void emitObjectRemoved(QObject *o);

    DeviceProfile deviceProfile() const;
    QString styleName() const;
    QString deviceProfileName() const;

    enum LineTerminatorMode {
        LFLineTerminator,
        CRLFLineTerminator,
        NativeLineTerminator =
#if defined (Q_OS_WIN)
            CRLFLineTerminator
#else
            LFLineTerminator
#endif
    };

    void setLineTerminatorMode(LineTerminatorMode mode);
    LineTerminatorMode lineTerminatorMode() const;

    bool useIdBasedTranslations() const;
    void setUseIdBasedTranslations(bool v);

    bool connectSlotsByName() const;
    void setConnectSlotsByName(bool v);

public slots:
    void resourceSetActivated(QtResourceSet *resourceSet, bool resourceSetChanged);

private slots:
    void triggerDefaultAction(QWidget *w);
    void sheetDestroyed(QObject *object);

private:
    void syncGridFeature();
    void connectSheet(QDesignerPropertySheet *sheet);
    void disconnectSheet(QDesignerPropertySheet *sheet);

    FormWindowBasePrivate *m_d;
};

}  // namespace qdesigner_internal

QT_END_NAMESPACE

#endif // FORMWINDOWBASE_H
