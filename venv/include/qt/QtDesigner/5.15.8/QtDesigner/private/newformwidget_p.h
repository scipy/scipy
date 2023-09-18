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

#ifndef NEWFORMWIDGET_H
#define NEWFORMWIDGET_H

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

#include "shared_global_p.h"
#include "deviceprofile_p.h"

#include <QtDesigner/abstractnewformwidget.h>

#include <QtWidgets/qwidget.h>
#include <QtGui/qpixmap.h>

#include <QtCore/qstringlist.h>
#include <QtCore/qpair.h>
#include <QtCore/qmap.h>
#include <QtCore/qvector.h>

QT_BEGIN_NAMESPACE

class QIODevice;
class QTreeWidgetItem;

namespace qdesigner_internal {

namespace Ui {
    class NewFormWidget;
}

class QDESIGNER_SHARED_EXPORT NewFormWidget : public QDesignerNewFormWidgetInterface
{
    Q_OBJECT
    Q_DISABLE_COPY_MOVE(NewFormWidget)

public:
    using DeviceProfileList = QVector<qdesigner_internal::DeviceProfile>;

    explicit NewFormWidget(QDesignerFormEditorInterface *core, QWidget *parentWidget);
    ~NewFormWidget() override;

    bool hasCurrentTemplate() const override;
    QString currentTemplate(QString *errorMessage = nullptr) override;

    // Convenience for implementing file dialogs with preview
    static QImage grabForm(QDesignerFormEditorInterface *core,
                           QIODevice &file,
                           const QString &workingDir,
                           const qdesigner_internal::DeviceProfile &dp);

private slots:
    void on_treeWidget_itemActivated(QTreeWidgetItem *item);
    void on_treeWidget_currentItemChanged(QTreeWidgetItem *current, QTreeWidgetItem *);
    void on_treeWidget_itemPressed(QTreeWidgetItem *item);
    void slotDeviceProfileIndexChanged(int idx);

private:
    QPixmap formPreviewPixmap(const QString &fileName) const;
    QPixmap formPreviewPixmap(QIODevice &file, const QString &workingDir = QString()) const;
    QPixmap formPreviewPixmap(const QTreeWidgetItem *item);

    void loadFrom(const QString &path, bool resourceFile, const QString &uiExtension,
                  const QString &selectedItem, QTreeWidgetItem *&selectedItemFound);
    void loadFrom(const QString &title, const QStringList &nameList,
                  const QString &selectedItem, QTreeWidgetItem *&selectedItemFound);

private:
    QString itemToTemplate(const QTreeWidgetItem *item, QString *errorMessage) const;
    QString currentTemplateI(QString *ptrToErrorMessage);

    QSize templateSize() const;
    void setTemplateSize(const QSize &s);
    int profileComboIndex() const;
    qdesigner_internal::DeviceProfile currentDeviceProfile() const;
    bool showCurrentItemPixmap();

    // Pixmap cache (item, profile combo index)
    using ItemPixmapCacheKey = QPair<const QTreeWidgetItem *, int>;
    using ItemPixmapCache = QMap<ItemPixmapCacheKey, QPixmap>;
    ItemPixmapCache m_itemPixmapCache;

    QDesignerFormEditorInterface *m_core;
    Ui::NewFormWidget *m_ui;
    QTreeWidgetItem *m_currentItem;
    QTreeWidgetItem *m_acceptedItem;
    DeviceProfileList m_deviceProfiles;
};

}

QT_END_NAMESPACE

#endif // NEWFORMWIDGET_H
