/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt Toolkit.
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

#ifndef QMEDIASERVICEPROVIDER_H
#define QMEDIASERVICEPROVIDER_H

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

#include <QtCore/qobject.h>
#include <QtCore/qshareddata.h>
#include <qtmultimediaglobal.h>
#include "qmultimedia.h"
#include "qmediaserviceproviderplugin.h"

QT_BEGIN_NAMESPACE


class QMediaService;

class Q_MULTIMEDIA_EXPORT QMediaServiceProvider : public QObject
{
    Q_OBJECT

public:
    virtual QMediaService* requestService(const QByteArray &type, const QMediaServiceProviderHint &hint = QMediaServiceProviderHint()) = 0;
    virtual void releaseService(QMediaService *service) = 0;

    virtual QMediaServiceProviderHint::Features supportedFeatures(const QMediaService *service) const;

    virtual QMultimedia::SupportEstimate hasSupport(const QByteArray &serviceType,
                                             const QString &mimeType,
                                             const QStringList& codecs,
                                             int flags = 0) const;
    virtual QStringList supportedMimeTypes(const QByteArray &serviceType, int flags = 0) const;

    virtual QByteArray defaultDevice(const QByteArray &serviceType) const;
    virtual QList<QByteArray> devices(const QByteArray &serviceType) const;
    virtual QString deviceDescription(const QByteArray &serviceType, const QByteArray &device);

    virtual QCamera::Position cameraPosition(const QByteArray &device) const;
    virtual int cameraOrientation(const QByteArray &device) const;

    static QMediaServiceProvider* defaultServiceProvider();
    static void setDefaultServiceProvider(QMediaServiceProvider *provider);
};

QT_END_NAMESPACE


#endif  // QMEDIASERVICEPROVIDER_H
