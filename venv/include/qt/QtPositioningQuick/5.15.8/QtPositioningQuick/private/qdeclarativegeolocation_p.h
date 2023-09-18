/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtPositioning module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:LGPL$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** GNU Lesser General Public License Usage
** Alternatively, this file may be used under the terms of the GNU Lesser
** General Public License version 3 as published by the Free Software
** Foundation and appearing in the file LICENSE.LGPL3 included in the
** packaging of this file. Please review the following information to
** ensure the GNU Lesser General Public License version 3 requirements
** will be met: https://www.gnu.org/licenses/lgpl-3.0.html.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 2.0 or (at your option) the GNU General
** Public license version 3 or any later version approved by the KDE Free
** Qt Foundation. The licenses are as published by the Free Software
** Foundation and appearing in the file LICENSE.GPL2 and LICENSE.GPL3
** included in the packaging of this file. Please review the following
** information to ensure the GNU General Public License requirements will
** be met: https://www.gnu.org/licenses/gpl-2.0.html and
** https://www.gnu.org/licenses/gpl-3.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef QDECLARATIVEGEOLOCATION_P_H
#define QDECLARATIVEGEOLOCATION_P_H

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

#include <QtCore/QObject>
#include <QtCore/QVariantMap>
#include <QtPositioning/QGeoLocation>
#include <QtPositioning/qgeorectangle.h>
#include <QtPositioningQuick/private/qdeclarativegeoaddress_p.h>
#include <QtPositioningQuick/private/qpositioningquickglobal_p.h>

QT_BEGIN_NAMESPACE

class Q_POSITIONINGQUICK_PRIVATE_EXPORT QDeclarativeGeoLocation : public QObject
{
    Q_OBJECT

    Q_PROPERTY(QGeoLocation location READ location WRITE setLocation)
    Q_PROPERTY(QDeclarativeGeoAddress *address READ address WRITE setAddress NOTIFY addressChanged)
    Q_PROPERTY(QGeoCoordinate coordinate READ coordinate WRITE setCoordinate NOTIFY coordinateChanged)
    Q_PROPERTY(QGeoRectangle boundingBox READ boundingBox WRITE setBoundingBox NOTIFY boundingBoxChanged)
    Q_PROPERTY(QVariantMap extendedAttributes MEMBER m_extendedAttributes NOTIFY extendedAttributesChanged REVISION 13)

public:
    explicit QDeclarativeGeoLocation(QObject *parent = 0);
    explicit QDeclarativeGeoLocation(const QGeoLocation &src, QObject *parent = 0);
    ~QDeclarativeGeoLocation();

    QGeoLocation location() const;
    void setLocation(const QGeoLocation &src);

    QDeclarativeGeoAddress *address() const;
    void setAddress(QDeclarativeGeoAddress *address);
    QGeoCoordinate coordinate() const;
    void setCoordinate(const QGeoCoordinate coordinate);

    QGeoRectangle boundingBox() const;
    void setBoundingBox(const QGeoRectangle &boundingBox);

Q_SIGNALS:
    void addressChanged();
    void coordinateChanged();
    void boundingBoxChanged();
    void extendedAttributesChanged();

private:
    QDeclarativeGeoAddress *m_address = nullptr;
    QGeoRectangle m_boundingBox;
    QGeoCoordinate m_coordinate;
    QVariantMap m_extendedAttributes;
};

QT_END_NAMESPACE

#endif // QDECLARATIVELOCATION_P_H
