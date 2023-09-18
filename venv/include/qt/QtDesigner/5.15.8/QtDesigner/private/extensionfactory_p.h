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

#ifndef SHARED_EXTENSIONFACTORY_H
#define SHARED_EXTENSIONFACTORY_H

#include <QtDesigner/default_extensionfactory.h>
#include <QtDesigner/qextensionmanager.h>

QT_BEGIN_NAMESPACE

namespace qdesigner_internal {

// Extension factory for registering an extension for an object type.
template <class ExtensionInterface, class Object, class Extension>
class ExtensionFactory: public QExtensionFactory
{
public:
    explicit ExtensionFactory(const QString &iid, QExtensionManager *parent = nullptr);

    // Convenience for registering the extension. Do not use for derived classes.
    static void registerExtension(QExtensionManager *mgr, const QString &iid);

protected:
    QObject *createExtension(QObject *qObject, const QString &iid, QObject *parent) const override;

private:
    // Can be overwritten to perform checks on the object.
    // Default does a qobject_cast to the desired class.
    virtual Object *checkObject(QObject *qObject) const;

    const QString m_iid;
};

template <class ExtensionInterface, class Object, class Extension>
ExtensionFactory<ExtensionInterface, Object, Extension>::ExtensionFactory(const QString &iid, QExtensionManager *parent) :
    QExtensionFactory(parent),
    m_iid(iid)
{
}

template <class ExtensionInterface, class Object, class Extension>
Object *ExtensionFactory<ExtensionInterface, Object, Extension>::checkObject(QObject *qObject) const
{
    return qobject_cast<Object*>(qObject);
}

template <class ExtensionInterface, class Object, class Extension>
QObject *ExtensionFactory<ExtensionInterface, Object, Extension>::createExtension(QObject *qObject, const QString &iid, QObject *parent) const
{
    if (iid != m_iid)
        return nullptr;

    Object *object = checkObject(qObject);
    if (!object)
        return nullptr;

    return new Extension(object, parent);
}

template <class ExtensionInterface, class Object, class Extension>
void ExtensionFactory<ExtensionInterface, Object, Extension>::registerExtension(QExtensionManager *mgr, const QString &iid)
{
    ExtensionFactory *factory = new ExtensionFactory(iid, mgr);
    mgr->registerExtensions(factory, iid);
}
}  // namespace qdesigner_internal

QT_END_NAMESPACE

#endif // SHARED_EXTENSIONFACTORY_H
