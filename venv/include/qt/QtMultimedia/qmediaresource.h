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

#ifndef QMEDIARESOURCE_H
#define QMEDIARESOURCE_H

#include <QtCore/qmap.h>
#include <QtCore/qmetatype.h>
#include <QtNetwork/qnetworkrequest.h>

#include <QtMultimedia/qtmultimediaglobal.h>

QT_BEGIN_NAMESPACE

// Class forward declaration required for QDoc bug
class QString;
class Q_MULTIMEDIA_EXPORT QMediaResource
{
public:
    QMediaResource();
    QMediaResource(const QUrl &url, const QString &mimeType = QString());
    QMediaResource(const QNetworkRequest &request, const QString &mimeType = QString());
    QMediaResource(const QMediaResource &other);
    QMediaResource &operator =(const QMediaResource &other);
    ~QMediaResource();

    bool isNull() const;

    bool operator ==(const QMediaResource &other) const;
    bool operator !=(const QMediaResource &other) const;

    QUrl url() const;
    QNetworkRequest request() const;
    QString mimeType() const;

    QString language() const;
    void setLanguage(const QString &language);

    QString audioCodec() const;
    void setAudioCodec(const QString &codec);

    QString videoCodec() const;
    void setVideoCodec(const QString &codec);

    qint64 dataSize() const;
    void setDataSize(const qint64 size);

    int audioBitRate() const;
    void setAudioBitRate(int rate);

    int sampleRate() const;
    void setSampleRate(int frequency);

    int channelCount() const;
    void setChannelCount(int channels);

    int videoBitRate() const;
    void setVideoBitRate(int rate);

    QSize resolution() const;
    void setResolution(const QSize &resolution);
    void setResolution(int width, int height);


private:
    enum Property
    {
        Url,
        Request,
        MimeType,
        Language,
        AudioCodec,
        VideoCodec,
        DataSize,
        AudioBitRate,
        VideoBitRate,
        SampleRate,
        ChannelCount,
        Resolution
    };
    QMap<int, QVariant> values;
};

typedef QList<QMediaResource> QMediaResourceList;

QT_END_NAMESPACE

Q_DECLARE_METATYPE(QMediaResource)
Q_DECLARE_METATYPE(QMediaResourceList)

#endif
