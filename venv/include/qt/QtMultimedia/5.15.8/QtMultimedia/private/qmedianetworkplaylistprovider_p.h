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

#ifndef QMEDIANETWORKPAYLISTPROVIDER_P_H
#define QMEDIANETWORKPAYLISTPROVIDER_P_H

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API. It exists purely as an
// implementation detail. This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.
//

#include "qmediaplaylistprovider_p.h"

QT_BEGIN_NAMESPACE


class QMediaNetworkPlaylistProviderPrivate;
class Q_MULTIMEDIA_EXPORT QMediaNetworkPlaylistProvider : public QMediaPlaylistProvider
{
    Q_OBJECT
public:
    QMediaNetworkPlaylistProvider(QObject *parent = nullptr);
    ~QMediaNetworkPlaylistProvider();

    bool load(const QNetworkRequest &request, const char *format = nullptr) override;

    int mediaCount() const override;
    QMediaContent media(int pos) const override;

    bool isReadOnly() const override;

    bool addMedia(const QMediaContent &content) override;
    bool addMedia(const QList<QMediaContent> &items) override;
    bool insertMedia(int pos, const QMediaContent &content) override;
    bool insertMedia(int pos, const QList<QMediaContent> &items) override;
    bool moveMedia(int from, int to) override;
    bool removeMedia(int pos) override;
    bool removeMedia(int start, int end) override;
    bool clear() override;

public Q_SLOTS:
    void shuffle() override;

private:
    Q_DISABLE_COPY(QMediaNetworkPlaylistProvider)
    Q_DECLARE_PRIVATE(QMediaNetworkPlaylistProvider)
    Q_PRIVATE_SLOT(d_func(), void _q_handleParserError(QPlaylistFileParser::ParserError err, const QString &))
    Q_PRIVATE_SLOT(d_func(), void _q_handleNewItem(const QVariant& content))
};

QT_END_NAMESPACE


#endif // QMEDIANETWORKPAYLISTSOURCE_P_H
