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

#ifndef QRADIOTUNERCONTROL_H
#define QRADIOTUNERCONTROL_H

#include <QtMultimedia/qmediacontrol.h>
#include <QtMultimedia/qradiotuner.h>

QT_BEGIN_NAMESPACE

// Required for QDoc workaround
class QString;

class Q_MULTIMEDIA_EXPORT QRadioTunerControl : public QMediaControl
{
    Q_OBJECT

public:
    ~QRadioTunerControl();

    virtual QRadioTuner::State state() const = 0;

    virtual QRadioTuner::Band band() const = 0;
    virtual void setBand(QRadioTuner::Band b) = 0;
    virtual bool isBandSupported(QRadioTuner::Band b) const = 0;

    virtual int frequency() const = 0;
    virtual int frequencyStep(QRadioTuner::Band b) const = 0;
    virtual QPair<int,int> frequencyRange(QRadioTuner::Band b) const = 0;
    virtual void setFrequency(int frequency) = 0;

    virtual bool isStereo() const = 0;
    virtual QRadioTuner::StereoMode stereoMode() const = 0;
    virtual void setStereoMode(QRadioTuner::StereoMode mode) = 0;

    virtual int signalStrength() const = 0;

    virtual int volume() const = 0;
    virtual void setVolume(int volume) = 0;

    virtual bool isMuted() const = 0;
    virtual void setMuted(bool muted) = 0;

    virtual bool isSearching() const = 0;

    virtual bool isAntennaConnected() const { return true; }

    virtual void searchForward() = 0;
    virtual void searchBackward() = 0;
    virtual void searchAllStations(QRadioTuner::SearchMode searchMode = QRadioTuner::SearchFast) = 0;
    virtual void cancelSearch() = 0;

    virtual void start() = 0;
    virtual void stop() = 0;

    virtual QRadioTuner::Error error() const = 0;
    virtual QString errorString() const = 0;

Q_SIGNALS:
    void stateChanged(QRadioTuner::State state);
    void bandChanged(QRadioTuner::Band band);
    void frequencyChanged(int frequency);
    void stereoStatusChanged(bool stereo);
    void searchingChanged(bool searching);
    void signalStrengthChanged(int signalStrength);
    void volumeChanged(int volume);
    void mutedChanged(bool muted);
    void error(QRadioTuner::Error err);
    void stationFound(int frequency, QString stationId);
    void antennaConnectedChanged(bool connectionStatus);

protected:
    explicit QRadioTunerControl(QObject *parent = nullptr);
};

#define QRadioTunerControl_iid "org.qt-project.qt.radiotunercontrol/5.0"
Q_MEDIA_DECLARE_CONTROL(QRadioTunerControl, QRadioTunerControl_iid)

QT_END_NAMESPACE


#endif  // QRADIOTUNERCONTROL_H
