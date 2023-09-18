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

#ifndef ABSTRACTPROMOTIONINTERFACE_H
#define ABSTRACTPROMOTIONINTERFACE_H

#include <QtDesigner/sdk_global.h>

#include <QtCore/qpair.h>
#include <QtCore/qlist.h>
#include <QtCore/qset.h>

QT_BEGIN_NAMESPACE

class QDesignerWidgetDataBaseItemInterface;

class QDESIGNER_SDK_EXPORT QDesignerPromotionInterface
{
public:
    virtual ~QDesignerPromotionInterface();

    struct PromotedClass {
        QDesignerWidgetDataBaseItemInterface *baseItem;
        QDesignerWidgetDataBaseItemInterface *promotedItem;
    };

    using PromotedClasses = QList<PromotedClass>;

    virtual PromotedClasses promotedClasses() const = 0;

    virtual QSet<QString> referencedPromotedClassNames()  const = 0;

    virtual bool addPromotedClass(const QString &baseClass,
                                  const QString &className,
                                  const QString &includeFile,
                                  QString *errorMessage) = 0;

    virtual bool removePromotedClass(const QString &className, QString *errorMessage) = 0;

    virtual bool changePromotedClassName(const QString &oldClassName, const QString &newClassName, QString *errorMessage) = 0;

    virtual bool setPromotedClassIncludeFile(const QString &className, const QString &includeFile, QString *errorMessage) = 0;

    virtual QList<QDesignerWidgetDataBaseItemInterface *> promotionBaseClasses() const = 0;
};

QT_END_NAMESPACE

#endif // ABSTRACTPROMOTIONINTERFACE_H
