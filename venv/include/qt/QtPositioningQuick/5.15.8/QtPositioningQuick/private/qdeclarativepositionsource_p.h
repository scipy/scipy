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
***************************************************************************/

#ifndef QDECLARATIVEPOSITIONSOURCE_H
#define QDECLARATIVEPOSITIONSOURCE_H

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

#include <QtPositioningQuick/private/qpositioningquickglobal_p.h>
#include <QtPositioningQuick/private/qdeclarativeposition_p.h>
#include <QtCore/QObject>
#include <QtNetwork/QAbstractSocket>
#include <QtQml/QQmlParserStatus>
#include <QtPositioning/qgeopositioninfosource.h>
#include <QtPositioningQuick/private/qdeclarativepluginparameter_p.h>

QT_BEGIN_NAMESPACE

class QFile;
class QTcpSocket;

class Q_POSITIONINGQUICK_PRIVATE_EXPORT QDeclarativePositionSource : public QObject, public QQmlParserStatus
{
    Q_OBJECT

    Q_PROPERTY(QDeclarativePosition *position READ position NOTIFY positionChanged)
    Q_PROPERTY(bool active READ isActive WRITE setActive NOTIFY activeChanged)
    Q_PROPERTY(bool valid READ isValid NOTIFY validityChanged)
    Q_PROPERTY(QUrl nmeaSource READ nmeaSource WRITE setNmeaSource NOTIFY nmeaSourceChanged)
    Q_PROPERTY(int updateInterval READ updateInterval WRITE setUpdateInterval NOTIFY updateIntervalChanged)
    Q_PROPERTY(PositioningMethods supportedPositioningMethods READ supportedPositioningMethods NOTIFY supportedPositioningMethodsChanged)
    Q_PROPERTY(PositioningMethods preferredPositioningMethods READ preferredPositioningMethods WRITE setPreferredPositioningMethods NOTIFY preferredPositioningMethodsChanged)
    Q_PROPERTY(SourceError sourceError READ sourceError NOTIFY sourceErrorChanged)
    Q_PROPERTY(QString name READ name WRITE setName NOTIFY nameChanged)
    Q_PROPERTY(QQmlListProperty<QDeclarativePluginParameter> parameters READ parameters REVISION 14)
    Q_ENUMS(PositioningMethod)

    Q_CLASSINFO("DefaultProperty", "parameters")
    Q_INTERFACES(QQmlParserStatus)

public:
    enum PositioningMethod {
        NoPositioningMethods = QGeoPositionInfoSource::NoPositioningMethods,
        SatellitePositioningMethods = QGeoPositionInfoSource::SatellitePositioningMethods,
        NonSatellitePositioningMethods = QGeoPositionInfoSource::NonSatellitePositioningMethods,
        AllPositioningMethods = QGeoPositionInfoSource::AllPositioningMethods
    };

    Q_DECLARE_FLAGS(PositioningMethods, PositioningMethod)
    Q_FLAGS(PositioningMethods)

    enum SourceError {
        AccessError = QGeoPositionInfoSource::AccessError,
        ClosedError = QGeoPositionInfoSource::ClosedError,
        UnknownSourceError = QGeoPositionInfoSource::UnknownSourceError,
        NoError = QGeoPositionInfoSource::NoError,

        //Leave a gap for future error enum values in QGeoPositionInfoSource::Error
        SocketError = 100
    };
    Q_ENUMS(SourceError)

    QDeclarativePositionSource();
    ~QDeclarativePositionSource();
    void setNmeaSource(const QUrl &nmeaSource);
    void setUpdateInterval(int updateInterval);
    void setActive(bool active);
    void setPreferredPositioningMethods(PositioningMethods methods);

    QString name() const;
    void setName(const QString &name);

    QUrl nmeaSource() const;
    int updateInterval() const;
    bool isActive() const;
    bool isValid() const;
    QDeclarativePosition *position();
    PositioningMethods supportedPositioningMethods() const;
    PositioningMethods preferredPositioningMethods() const;
    SourceError sourceError() const;
    QGeoPositionInfoSource *positionSource() const;
    QQmlListProperty<QDeclarativePluginParameter> parameters();
    QVariantMap parameterMap() const;

    // Virtuals from QQmlParserStatus
    void classBegin() { }
    void componentComplete();

    Q_REVISION(14) Q_INVOKABLE bool setBackendProperty(const QString &name, const QVariant &value);
    Q_REVISION(14) Q_INVOKABLE QVariant backendProperty(const QString &name) const;

public Q_SLOTS:
    void update(); // TODO Qt 6 change to void update(int)
    void start();
    void stop();

Q_SIGNALS:
    void positionChanged();
    void activeChanged();
    void nmeaSourceChanged();
    void updateIntervalChanged();
    void supportedPositioningMethodsChanged();
    void preferredPositioningMethodsChanged();
    void sourceErrorChanged();
    void nameChanged();
    void validityChanged();
    void updateTimeout();

private Q_SLOTS:
    void positionUpdateReceived(const QGeoPositionInfo &update);
    void sourceErrorReceived(const QGeoPositionInfoSource::Error error);
    void socketConnected();
    void socketError(QAbstractSocket::SocketError error);
    void updateTimeoutReceived();
    void onParameterInitialized();

private:
    void setPosition(const QGeoPositionInfo &pi);
    void setSource(QGeoPositionInfoSource *source);
    bool parametersReady();
    void tryAttach(const QString &name, bool useFallback = true);

    static void parameter_append(QQmlListProperty<QDeclarativePluginParameter> *prop, QDeclarativePluginParameter *mapObject);
    static int parameter_count(QQmlListProperty<QDeclarativePluginParameter> *prop);
    static QDeclarativePluginParameter *parameter_at(QQmlListProperty<QDeclarativePluginParameter> *prop, int index);
    static void parameter_clear(QQmlListProperty<QDeclarativePluginParameter> *prop);

    QGeoPositionInfoSource *m_positionSource;
    QDeclarativePosition m_position;
    PositioningMethods m_preferredPositioningMethods;
    QFile *m_nmeaFile;
    QTcpSocket *m_nmeaSocket;
    QString m_nmeaFileName;
    QUrl m_nmeaSource;
    QString m_providerName;
    bool m_active;
    bool m_singleUpdate;
    int m_updateInterval;
    SourceError m_sourceError;
    QList<QDeclarativePluginParameter *> m_parameters;
    bool m_componentComplete = false;
    bool m_parametersInitialized = false;
};

QT_END_NAMESPACE

QML_DECLARE_TYPE(QDeclarativePositionSource)

#endif
