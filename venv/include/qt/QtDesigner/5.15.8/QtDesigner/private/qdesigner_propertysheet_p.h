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

#ifndef QDESIGNER_PROPERTYSHEET_H
#define QDESIGNER_PROPERTYSHEET_H

#include "shared_global_p.h"
#include "dynamicpropertysheet.h"
#include <QtDesigner/propertysheet.h>
#include <QtDesigner/default_extensionfactory.h>
#include <QtDesigner/qextensionmanager.h>

#include <QtCore/qvariant.h>
#include <QtCore/qpair.h>

#include <QtCore/qpointer.h>

QT_BEGIN_NAMESPACE

class QLayout;
class QDesignerFormEditorInterface;
class QDesignerPropertySheetPrivate;

namespace qdesigner_internal
{
    class DesignerPixmapCache;
    class DesignerIconCache;
    class FormWindowBase;
}

class QDESIGNER_SHARED_EXPORT QDesignerPropertySheet: public QObject, public QDesignerPropertySheetExtension, public QDesignerDynamicPropertySheetExtension
{
    Q_OBJECT
    Q_INTERFACES(QDesignerPropertySheetExtension QDesignerDynamicPropertySheetExtension)
public:
    explicit QDesignerPropertySheet(QObject *object, QObject *parent = nullptr);
    ~QDesignerPropertySheet() override;

    int indexOf(const QString &name) const override;

    int count() const override;
    QString propertyName(int index) const override;

    QString propertyGroup(int index) const override;
    void setPropertyGroup(int index, const QString &group) override;

    bool hasReset(int index) const override;
    bool reset(int index) override;

    bool isAttribute(int index) const override;
    void setAttribute(int index, bool b) override;

    bool isVisible(int index) const override;
    void setVisible(int index, bool b) override;

    QVariant property(int index) const override;
    void setProperty(int index, const QVariant &value) override;

    bool isChanged(int index) const override;

    void setChanged(int index, bool changed) override;

    bool dynamicPropertiesAllowed() const override;
    int addDynamicProperty(const QString &propertyName, const QVariant &value) override;
    bool removeDynamicProperty(int index) override;
    bool isDynamicProperty(int index) const override;
    bool canAddDynamicProperty(const QString &propertyName) const override;

    bool isDefaultDynamicProperty(int index) const;

    bool isResourceProperty(int index) const;
    QVariant defaultResourceProperty(int index) const;

    qdesigner_internal::DesignerPixmapCache *pixmapCache() const;
    void setPixmapCache(qdesigner_internal::DesignerPixmapCache *cache);
    qdesigner_internal::DesignerIconCache *iconCache() const;
    void setIconCache(qdesigner_internal::DesignerIconCache *cache);
    int createFakeProperty(const QString &propertyName, const QVariant &value = QVariant());

    bool isEnabled(int index) const override;
    QObject *object() const;

    static bool internalDynamicPropertiesEnabled();
    static void setInternalDynamicPropertiesEnabled(bool v);

protected:
    bool isAdditionalProperty(int index) const;
    bool isFakeProperty(int index) const;
    QVariant resolvePropertyValue(int index, const QVariant &value) const;
    QVariant metaProperty(int index) const;
    void setFakeProperty(int index, const QVariant &value);
    void clearFakeProperties();

    bool isFakeLayoutProperty(int index) const;
    bool isDynamic(int index) const;
    qdesigner_internal::FormWindowBase *formWindowBase() const;
    QDesignerFormEditorInterface *core() const;

public:
    enum PropertyType { PropertyNone,
                        PropertyLayoutObjectName,
                        PropertyLayoutLeftMargin,
                        PropertyLayoutTopMargin,
                        PropertyLayoutRightMargin,
                        PropertyLayoutBottomMargin,
                        PropertyLayoutSpacing,
                        PropertyLayoutHorizontalSpacing,
                        PropertyLayoutVerticalSpacing,
                        PropertyLayoutSizeConstraint,
                        PropertyLayoutFieldGrowthPolicy,
                        PropertyLayoutRowWrapPolicy,
                        PropertyLayoutLabelAlignment,
                        PropertyLayoutFormAlignment,
                        PropertyLayoutBoxStretch,
                        PropertyLayoutGridRowStretch,
                        PropertyLayoutGridColumnStretch,
                        PropertyLayoutGridRowMinimumHeight,
                        PropertyLayoutGridColumnMinimumWidth,
                        PropertyBuddy,
                        PropertyAccessibility,
                        PropertyGeometry,
                        PropertyChecked,
                        PropertyCheckable,
                        PropertyVisible,
                        PropertyWindowTitle,
                        PropertyWindowIcon,
                        PropertyWindowFilePath,
                        PropertyWindowOpacity,
                        PropertyWindowIconText,
                        PropertyWindowModality,
                        PropertyWindowModified,
                        PropertyStyleSheet,
                        PropertyText
    };

    enum ObjectType { ObjectNone, ObjectLabel, ObjectLayout, ObjectLayoutWidget };
    enum ObjectFlag
    {
        CheckableProperty = 0x1 // Has a "checked" property depending on "checkable"
    };
    Q_DECLARE_FLAGS(ObjectFlags, ObjectFlag)

    static ObjectType objectTypeFromObject(const QObject *o);
    static ObjectFlags objectFlagsFromObject(const QObject *o);
    static PropertyType propertyTypeFromName(const QString &name);

protected:
    PropertyType propertyType(int index) const;
    ObjectType objectType() const;

private:
    QDesignerPropertySheetPrivate *d;
};

/* Abstract base class for factories that register a property sheet that implements
 * both QDesignerPropertySheetExtension and QDesignerDynamicPropertySheetExtension
 * by multiple inheritance. The factory maintains ownership of
 * the extension and returns it for both id's. */

class QDESIGNER_SHARED_EXPORT QDesignerAbstractPropertySheetFactory: public QExtensionFactory
{
    Q_OBJECT
    Q_INTERFACES(QAbstractExtensionFactory)
public:
    explicit QDesignerAbstractPropertySheetFactory(QExtensionManager *parent = nullptr);
    ~QDesignerAbstractPropertySheetFactory() override;

    QObject *extension(QObject *object, const QString &iid) const override;

private slots:
    void objectDestroyed(QObject *object);

private:
    virtual QObject *createPropertySheet(QObject *qObject, QObject *parent) const = 0;

    struct PropertySheetFactoryPrivate;
    PropertySheetFactoryPrivate *m_impl;
};

/* Convenience factory template for property sheets that implement
 * QDesignerPropertySheetExtension and QDesignerDynamicPropertySheetExtension
 * by multiple inheritance. */

template <class Object, class PropertySheet>
class QDesignerPropertySheetFactory : public QDesignerAbstractPropertySheetFactory {
public:
    explicit QDesignerPropertySheetFactory(QExtensionManager *parent = nullptr);

    static void registerExtension(QExtensionManager *mgr);

private:
    // Does a  qobject_cast on  the object.
    QObject *createPropertySheet(QObject *qObject, QObject *parent) const override;
};

template <class Object, class PropertySheet>
QDesignerPropertySheetFactory<Object, PropertySheet>::QDesignerPropertySheetFactory(QExtensionManager *parent) :
    QDesignerAbstractPropertySheetFactory(parent)
{
}

template <class Object, class PropertySheet>
QObject *QDesignerPropertySheetFactory<Object, PropertySheet>::createPropertySheet(QObject *qObject, QObject *parent) const
{
    Object *object = qobject_cast<Object *>(qObject);
    if (!object)
        return nullptr;
    return new PropertySheet(object, parent);
}

template <class Object, class PropertySheet>
void QDesignerPropertySheetFactory<Object, PropertySheet>::registerExtension(QExtensionManager *mgr)
{
    QDesignerPropertySheetFactory *factory = new QDesignerPropertySheetFactory(mgr);
    mgr->registerExtensions(factory, Q_TYPEID(QDesignerPropertySheetExtension));
    mgr->registerExtensions(factory, Q_TYPEID(QDesignerDynamicPropertySheetExtension));
}


// Standard property sheet
typedef QDesignerPropertySheetFactory<QObject, QDesignerPropertySheet> QDesignerDefaultPropertySheetFactory;

Q_DECLARE_OPERATORS_FOR_FLAGS(QDesignerPropertySheet::ObjectFlags)

QT_END_NAMESPACE

#endif // QDESIGNER_PROPERTYSHEET_H
