/****************************************************************************
**
** Copyright (C) 2017 Klaralvdalens Datakonsult AB (KDAB).
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt3D module of the Qt Toolkit.
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

#ifndef QT3DCORE_QDOWNLOADNETWORKWORKER_P_H
#define QT3DCORE_QDOWNLOADNETWORKWORKER_P_H

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API.  It exists for the convenience
// of other Qt classes.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.
//

#include <QSharedPointer>
#include <QObject>
#include <QUrl>
#include <Qt3DCore/private/qdownloadhelperservice_p.h>
#include <QMutex>

QT_BEGIN_NAMESPACE

class QThread;
class QNetworkAccessManager;
class QNetworkReply;

namespace Qt3DCore {

class QDownloadRequest;
typedef QSharedPointer<QDownloadRequest> QDownloadRequestPtr;

class QDownloadNetworkWorker : public QObject
{
    Q_OBJECT
public:
    QDownloadNetworkWorker(QObject *parent = nullptr);

signals:
    void submitRequest(const Qt3DCore::QDownloadRequestPtr &request);
    void cancelRequest(const Qt3DCore::QDownloadRequestPtr &request);
    void cancelAllRequests();

    void requestDownloaded(const Qt3DCore::QDownloadRequestPtr &request);

private Q_SLOTS:
    void onRequestSubmited(const Qt3DCore::QDownloadRequestPtr &request);
    void onRequestCancelled(const Qt3DCore::QDownloadRequestPtr &request);
    void onAllRequestsCancelled();
    void onRequestFinished(QNetworkReply *reply);
    void onDownloadProgressed(qint64 bytesReceived, qint64 bytesTotal);

private:
    QNetworkAccessManager *m_networkManager;
    QVector< QPair<QDownloadRequestPtr, QNetworkReply *> > m_requests;
    QMutex m_mutex;
};

} // namespace Qt3DCore

QT_END_NAMESPACE

#endif // QT3DCORE_QDOWNLOADNETWORKWORKER_P_H
