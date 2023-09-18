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

#ifndef QRADIOTUNER_H
#define QRADIOTUNER_H

#include <QtCore/qobject.h>

#include <QtMultimedia/qmediaobject.h>
#include <QtMultimedia/qmediaenumdebug.h>

#include <QtCore/qpair.h>

QT_BEGIN_NAMESPACE

class QRadioData;
class QRadioTunerPrivate;
class Q_MULTIMEDIA_EXPORT QRadioTuner : public QMediaObject
{
    Q_OBJECT
    Q_PROPERTY(State state READ state NOTIFY stateChanged)
    Q_PROPERTY(Band band READ band WRITE setBand NOTIFY bandChanged)
    Q_PROPERTY(int frequency READ frequency WRITE setFrequency NOTIFY frequencyChanged)
    Q_PROPERTY(bool stereo READ isStereo NOTIFY stereoStatusChanged)
    Q_PROPERTY(StereoMode stereoMode READ stereoMode WRITE setStereoMode)
    Q_PROPERTY(int signalStrength READ signalStrength NOTIFY signalStrengthChanged)
    Q_PROPERTY(int volume READ volume WRITE setVolume NOTIFY volumeChanged)
    Q_PROPERTY(bool muted READ isMuted WRITE setMuted NOTIFY mutedChanged)
    Q_PROPERTY(bool searching READ isSearching NOTIFY searchingChanged)
    Q_PROPERTY(bool antennaConnected READ isAntennaConnected NOTIFY antennaConnectedChanged)
    Q_PROPERTY(QRadioData *radioData READ radioData CONSTANT)
    Q_ENUMS(State)
    Q_ENUMS(Band)
    Q_ENUMS(Error)
    Q_ENUMS(StereoMode)
    Q_ENUMS(SearchMode)

public:
    enum State { ActiveState, StoppedState };
    enum Band { AM, FM, SW, LW, FM2 };
    enum Error { NoError, ResourceError, OpenError, OutOfRangeError };
    enum StereoMode { ForceStereo, ForceMono, Auto };
    enum SearchMode { SearchFast, SearchGetStationId };

    explicit QRadioTuner(QObject *parent = nullptr);
    ~QRadioTuner();

    QMultimedia::AvailabilityStatus availability() const override;

    State state() const;

    Band band() const;

    bool isBandSupported(Band b) const;

    int frequency() const;
    int frequencyStep(Band band) const;
    QPair<int,int> frequencyRange(Band band) const;

    bool isStereo() const;
    void setStereoMode(QRadioTuner::StereoMode mode);
    StereoMode stereoMode() const;

    int signalStrength() const;

    int volume() const;
    bool isMuted() const;

    bool isSearching() const;

    bool isAntennaConnected() const;

    Error error() const;
    QString errorString() const;

    QRadioData *radioData() const;

public Q_SLOTS:
    void searchForward();
    void searchBackward();
    void searchAllStations(QRadioTuner::SearchMode searchMode = QRadioTuner::SearchFast);
    void cancelSearch();

    void setBand(Band band);
    void setFrequency(int frequency);

    void setVolume(int volume);
    void setMuted(bool muted);

    void start();
    void stop();

Q_SIGNALS:
    void stateChanged(QRadioTuner::State state);
    void bandChanged(QRadioTuner::Band band);
    void frequencyChanged(int frequency);
    void stereoStatusChanged(bool stereo);
    void searchingChanged(bool searching);
    void signalStrengthChanged(int signalStrength);
    void volumeChanged(int volume);
    void mutedChanged(bool muted);
    void stationFound(int frequency, QString stationId);
    void antennaConnectedChanged(bool connectionStatus);

    void error(QRadioTuner::Error error);

private:
    Q_DISABLE_COPY(QRadioTuner)
    Q_DECLARE_PRIVATE(QRadioTuner)
};

QT_END_NAMESPACE

Q_DECLARE_METATYPE(QRadioTuner::State)
Q_DECLARE_METATYPE(QRadioTuner::Band)
Q_DECLARE_METATYPE(QRadioTuner::Error)
Q_DECLARE_METATYPE(QRadioTuner::StereoMode)
Q_DECLARE_METATYPE(QRadioTuner::SearchMode)

Q_MEDIA_ENUM_DEBUG(QRadioTuner, State)
Q_MEDIA_ENUM_DEBUG(QRadioTuner, Band)
Q_MEDIA_ENUM_DEBUG(QRadioTuner, Error)
Q_MEDIA_ENUM_DEBUG(QRadioTuner, StereoMode)
Q_MEDIA_ENUM_DEBUG(QRadioTuner, SearchMode)

#endif  // QRADIOPLAYER_H
