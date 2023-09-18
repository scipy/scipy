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


#ifndef DESIGNERPROPERTYEDITOR_H
#define DESIGNERPROPERTYEDITOR_H

#include "shared_global_p.h"
#include "shared_enums_p.h"
#include <QtDesigner/abstractpropertyeditor.h>

QT_BEGIN_NAMESPACE

namespace qdesigner_internal {

// Extends the QDesignerPropertyEditorInterface by property comment handling and
// a signal for resetproperty.

class QDESIGNER_SHARED_EXPORT QDesignerPropertyEditor: public QDesignerPropertyEditorInterface
{
    Q_OBJECT
public:
    explicit QDesignerPropertyEditor(QWidget *parent = nullptr, Qt::WindowFlags flags = {});

    // A pair <ValidationMode, bool isTranslatable>.
    using StringPropertyParameters = QPair<TextPropertyValidationMode, bool>;

    // Return a pair of validation mode and flag indicating whether property is translatable
    // for textual properties.
    static StringPropertyParameters textPropertyValidationMode(QDesignerFormEditorInterface *core,
                const QObject *object, const QString &propertyName, bool isMainContainer);

Q_SIGNALS:
    void propertyValueChanged(const QString &name, const QVariant &value, bool enableSubPropertyHandling);
    void resetProperty(const QString &name);
    void addDynamicProperty(const QString &name, const QVariant &value);
    void removeDynamicProperty(const QString &name);
    void editorOpened();
    void editorClosed();

public Q_SLOTS:
    /* Quick update that assumes the actual count of properties has not changed
     * (as opposed to setObject()). N/A when for example executing a
     * layout command and margin properties appear. */
    virtual void updatePropertySheet() = 0;
    virtual void reloadResourceProperties() = 0;

private Q_SLOTS:
    void slotPropertyChanged(const QString &name, const QVariant &value);

protected:
    void emitPropertyValueChanged(const QString &name, const QVariant &value, bool enableSubPropertyHandling);

private:
    bool m_propertyChangedForwardingBlocked = false;

};

}  // namespace qdesigner_internal

QT_END_NAMESPACE

#endif // DESIGNERPROPERTYEDITOR_H
