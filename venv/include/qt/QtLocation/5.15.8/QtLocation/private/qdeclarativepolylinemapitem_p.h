/****************************************************************************
**
** Copyright (C) 2015 The Qt Company Ltd.
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

#ifndef QDECLARATIVEPOLYLINEMAPITEM
#define QDECLARATIVEPOLYLINEMAPITEM

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

#include <QtLocation/private/qlocationglobal_p.h>
#include <QtLocation/private/qdeclarativegeomapitembase_p.h>
#include <QtLocation/private/qgeomapitemgeometry_p.h>

#include <QtPositioning/QGeoPath>
#include <QtPositioning/private/qdoublevector2d_p.h>
#include <QSGGeometryNode>
#include <QSGFlatColorMaterial>

QT_BEGIN_NAMESPACE

class Q_LOCATION_PRIVATE_EXPORT QDeclarativeMapLineProperties : public QObject
{
    Q_OBJECT

    Q_PROPERTY(qreal width READ width WRITE setWidth NOTIFY widthChanged)
    Q_PROPERTY(QColor color READ color WRITE setColor NOTIFY colorChanged)

public:
    explicit QDeclarativeMapLineProperties(QObject *parent = 0);

    QColor color() const;
    void setColor(const QColor &color);

    qreal width() const;
    void setWidth(qreal width);

Q_SIGNALS:
    void widthChanged(qreal width);
    void colorChanged(const QColor &color);

private:
    qreal width_;
    QColor color_;
};

class QDeclarativePolylineMapItemPrivate;
class Q_LOCATION_PRIVATE_EXPORT QDeclarativePolylineMapItem : public QDeclarativeGeoMapItemBase
{
    Q_OBJECT
    Q_ENUMS(Backend)

    Q_PROPERTY(QJSValue path READ path WRITE setPath NOTIFY pathChanged)
    Q_PROPERTY(QDeclarativeMapLineProperties *line READ line CONSTANT)
    Q_PROPERTY(Backend backend READ backend WRITE setBackend NOTIFY backendChanged REVISION 15)

public:
    enum Backend {
        Software = 0,
#if QT_CONFIG(opengl)
        OpenGLLineStrip = 1,
        OpenGLExtruded = 2,
#endif
    };

    explicit QDeclarativePolylineMapItem(QQuickItem *parent = 0);
    ~QDeclarativePolylineMapItem();

    virtual void setMap(QDeclarativeGeoMap *quickMap, QGeoMap *map) override;
       //from QuickItem
    virtual QSGNode *updateMapItemPaintNode(QSGNode *, UpdatePaintNodeData *) override;

    Q_INVOKABLE int pathLength() const;
    Q_INVOKABLE void addCoordinate(const QGeoCoordinate &coordinate);
    Q_INVOKABLE void insertCoordinate(int index, const QGeoCoordinate &coordinate);
    Q_INVOKABLE void replaceCoordinate(int index, const QGeoCoordinate &coordinate);
    Q_INVOKABLE QGeoCoordinate coordinateAt(int index) const;
    Q_INVOKABLE bool containsCoordinate(const QGeoCoordinate &coordinate);
    Q_INVOKABLE void removeCoordinate(const QGeoCoordinate &coordinate);
    Q_INVOKABLE void removeCoordinate(int index);

    QJSValue path() const;
    virtual void setPath(const QJSValue &value);
    Q_INVOKABLE void setPath(const QGeoPath &path);

    bool contains(const QPointF &point) const override;
    const QGeoShape &geoShape() const override;
    void setGeoShape(const QGeoShape &shape) override;

    QDeclarativeMapLineProperties *line();

    Backend backend() const;
    void setBackend(Backend b);

Q_SIGNALS:
    void pathChanged();
    void backendChanged();

protected Q_SLOTS:
    void markSourceDirtyAndUpdate();
    void updateAfterLinePropertiesChanged();
    virtual void afterViewportChanged(const QGeoMapViewportChangeEvent &event) override;

protected:
    void geometryChanged(const QRectF &newGeometry, const QRectF &oldGeometry) override;
    void setPathFromGeoList(const QList<QGeoCoordinate> &path);
    void updatePolish() override;
    void componentComplete() override;
    void updateLineStyleParameter(QGeoMapParameter *p, const char *propertyName);
    void updateLineStyleParameter(QGeoMapParameter *p, const char *propertyName, bool update);

#ifdef QT_LOCATION_DEBUG
public:
#endif
    QGeoPath m_geopath;
    QDeclarativeMapLineProperties m_line;

    Backend m_backend = Software;
    bool m_dirtyMaterial;
    bool m_updatingGeometry;

    QScopedPointer<QDeclarativePolylineMapItemPrivate> m_d;

    friend class QDeclarativePolylineMapItemPrivate;
    friend class QDeclarativePolylineMapItemPrivateCPU;
    friend class QDeclarativePolylineMapItemPrivateOpenGLLineStrip;
    friend class QDeclarativePolylineMapItemPrivateOpenGLExtruded;
};

QT_END_NAMESPACE

QML_DECLARE_TYPE(QDeclarativeMapLineProperties)
QML_DECLARE_TYPE(QDeclarativePolylineMapItem)

#endif /* QDECLARATIVEPOLYLINEMAPITEM_H_ */
