/****************************************************************************
**
** Copyright (C) 2015 The Qt Company Ltd.
** Contact: http://www.qt.io/licensing/
**
** This file is part of the Purchasing module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:LGPL3-COMM$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see http://www.qt.io/terms-conditions. For further
** information use the contact form at http://www.qt.io/contact-us.
**
** GNU Lesser General Public License Usage
** Alternatively, this file may be used under the terms of the GNU Lesser
** General Public License version 3 as published by the Free Software
** Foundation and appearing in the file LICENSE.LGPLv3 included in the
** packaging of this file. Please review the following information to
** ensure the GNU Lesser General Public License version 3 requirements
** will be met: https://www.gnu.org/licenses/lgpl.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef QINAPPPRODUCT_H
#define QINAPPPRODUCT_H

#include <QtCore/qobject.h>
#include <QtCore/qsharedpointer.h>
#include <QtPurchasing/qtpurchasingglobal.h>

QT_BEGIN_NAMESPACE

class QInAppProductPrivate;
class Q_PURCHASING_EXPORT QInAppProduct: public QObject
{
    Q_OBJECT
    Q_ENUMS(ProductType)
    Q_PROPERTY(QString identifier READ identifier CONSTANT)
    Q_PROPERTY(ProductType productType READ productType CONSTANT)
    Q_PROPERTY(QString price READ price CONSTANT)
    Q_PROPERTY(QString title READ title CONSTANT)
    Q_PROPERTY(QString description READ description CONSTANT)

public:
    enum ProductType
    {
        Consumable,
        Unlockable
    };

    ~QInAppProduct();

    QString identifier() const;
    ProductType productType() const;

    QString price() const;
    QString title() const;
    QString description() const;

    Q_INVOKABLE virtual void purchase() = 0;

protected:
    explicit QInAppProduct(const QString &price, const QString &title, const QString &description, ProductType productType, const QString &identifier, QObject *parent = nullptr);

private:
    friend class QInAppStore;
    Q_DISABLE_COPY(QInAppProduct)

    QSharedPointer<QInAppProductPrivate> d;
};

QT_END_NAMESPACE

#endif // QINAPPPRODUCT_H
