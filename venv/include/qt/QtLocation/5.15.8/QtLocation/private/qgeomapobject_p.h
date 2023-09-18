/****************************************************************************
**
** Copyright (C) 2018 The Qt Company Ltd.
** Contact: http://www.qt.io/licensing/
**
** This file is part of the QtLocation module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:LGPL3$
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
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 2.0 or later as published by the Free
** Software Foundation and appearing in the file LICENSE.GPL included in
** the packaging of this file. Please review the following information to
** ensure the GNU General Public License version 2.0 requirements will be
** met: http://www.gnu.org/licenses/gpl-2.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef QGEOMAPOBJECTBASE_H
#define QGEOMAPOBJECTBASE_H

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

#include <QtLocation/private/qparameterizableobject_p.h>
#include <QExplicitlySharedDataPointer>
#include <QtPositioning/qgeoshape.h>
#include <qqml.h>

QT_BEGIN_NAMESPACE

class QGeoMapObjectPrivate;
class QGeoMap;

class Q_LOCATION_PRIVATE_EXPORT QGeoMapObject : public QParameterizableObject, public QQmlParserStatus
{
    Q_OBJECT

    Q_PROPERTY(bool visible READ visible WRITE setVisible NOTIFY visibleChanged)
    Q_PROPERTY(Type type READ type CONSTANT)
    Q_PROPERTY(QGeoShape geoShape READ geoShape WRITE setGeoShape STORED false) // non-NOTIFYable
    Q_INTERFACES(QQmlParserStatus)

public:
    enum Feature {
        NoFeature = 0x0,
        Clickable = 0x01,
        Draggable = 0x02,
        AllFeatures = 0xFFFFFFFF
    };

    enum Type {
        InvalidType = 0,
        ViewType = 1,
        RouteType = 2,
        RectangleType = 3,
        CircleType = 4,
        PolylineType = 5,
        PolygonType = 6,
        IconType = 7,
        UserType = 0x0100
    };

    Q_ENUM(Type)
    Q_DECLARE_FLAGS(Features, Feature)

    virtual ~QGeoMapObject();

    bool operator == (const QGeoMapObject &other) const;
    bool operator != (const QGeoMapObject &other) const;

    Features features() const;
    QGeoMapObjectPrivate *implementation() const;
    bool setImplementation(const QExplicitlySharedDataPointer<QGeoMapObjectPrivate> &pimpl);
    bool implemented() const;

    bool visible() const;
    void setVisible(bool visible);
    void setParentVisiblity(bool visible);

    Type type() const;

    virtual QList<QGeoMapObject*> geoMapObjectChildren() const;
    virtual void setMap(QGeoMap *map);
    QGeoMap *map() const;

    QGeoShape geoShape() const;
    void setGeoShape(const QGeoShape &shape);

Q_SIGNALS:
    void visibleChanged();
    void selected();
    void completed();

protected:
    QGeoMapObject(const QExplicitlySharedDataPointer<QGeoMapObjectPrivate> &dd, QObject *parent = nullptr);
    QExplicitlySharedDataPointer<QGeoMapObjectPrivate> d_ptr;

    void setChildrenVisibility();

    // QQmlParserStatus interface
    void classBegin() override;
    void componentComplete() override;
    void completeComponent();

    friend class QGeoMap;
    friend class QDeclarativeGeoMap;
    friend class QGeoMapLayer;
    friend class QDeclarativeNavigator;
};

QT_END_NAMESPACE

QML_DECLARE_TYPE(QGeoMapObject)

#endif // QGEOMAPOBJECTBASE_H
