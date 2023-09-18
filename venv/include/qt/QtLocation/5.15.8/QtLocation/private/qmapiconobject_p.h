/****************************************************************************
**
** Copyright (C) 2017 The Qt Company Ltd.
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

#ifndef QGEOMAPICON_P_H
#define QGEOMAPICON_P_H

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
#include <QtLocation/private/qgeomapobject_p.h>
#include <QtCore/QUrl>
#include <QGeoCoordinate>
#include <QtCore/qsize.h>

QT_BEGIN_NAMESPACE

class Q_LOCATION_PRIVATE_EXPORT QMapIconObject : public QGeoMapObject
{
    Q_OBJECT
    Q_PROPERTY(QGeoCoordinate coordinate READ coordinate WRITE setCoordinate NOTIFY coordinateChanged)
    Q_PROPERTY(QVariant content READ content WRITE setContent NOTIFY contentChanged)
    Q_PROPERTY(QSizeF iconSize READ iconSize WRITE setIconSize NOTIFY iconSizeChanged)

public:
    QMapIconObject(QObject *parent = nullptr);
    ~QMapIconObject() override;

    QVariant content() const;
    QGeoCoordinate coordinate() const;
    QSizeF iconSize() const;

    void setContent(QVariant content);
    void setCoordinate(const QGeoCoordinate &coordinate);
    void setIconSize(const QSizeF &size);

    void setMap(QGeoMap *map) override;

signals:
    void contentChanged(QVariant content);
    void coordinateChanged(QGeoCoordinate coordinate);
    void iconSizeChanged();
};

QT_END_NAMESPACE

#endif // QGEOMAPICON_P_H
