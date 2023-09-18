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

#ifndef QCAMERAIMAGEPROCESSING_H
#define QCAMERAIMAGEPROCESSING_H

#include <QtCore/qstringlist.h>
#include <QtCore/qpair.h>
#include <QtCore/qsize.h>
#include <QtCore/qpoint.h>
#include <QtCore/qrect.h>

#include <QtMultimedia/qmediacontrol.h>
#include <QtMultimedia/qmediaobject.h>
#include <QtMultimedia/qmediaservice.h>
#include <QtMultimedia/qmediaenumdebug.h>

QT_BEGIN_NAMESPACE


class QCamera;

class QCameraImageProcessingPrivate;
class Q_MULTIMEDIA_EXPORT QCameraImageProcessing : public QObject
{
    Q_OBJECT
    Q_ENUMS(WhiteBalanceMode ColorFilter)
public:
    enum WhiteBalanceMode {
        WhiteBalanceAuto = 0,
        WhiteBalanceManual = 1,
        WhiteBalanceSunlight = 2,
        WhiteBalanceCloudy = 3,
        WhiteBalanceShade = 4,
        WhiteBalanceTungsten = 5,
        WhiteBalanceFluorescent = 6,
        WhiteBalanceFlash = 7,
        WhiteBalanceSunset = 8,
        WhiteBalanceVendor = 1000
    };

    enum ColorFilter {
        ColorFilterNone,
        ColorFilterGrayscale,
        ColorFilterNegative,
        ColorFilterSolarize,
        ColorFilterSepia,
        ColorFilterPosterize,
        ColorFilterWhiteboard,
        ColorFilterBlackboard,
        ColorFilterAqua,
        ColorFilterVendor = 1000
    };

    bool isAvailable() const;

    WhiteBalanceMode whiteBalanceMode() const;
    void setWhiteBalanceMode(WhiteBalanceMode mode);
    bool isWhiteBalanceModeSupported(WhiteBalanceMode mode) const;

    qreal manualWhiteBalance() const;
    void setManualWhiteBalance(qreal colorTemperature);

    qreal brightness() const;
    void setBrightness(qreal value);

    qreal contrast() const;
    void setContrast(qreal value);

    qreal saturation() const;
    void setSaturation(qreal value);

    qreal sharpeningLevel() const;
    void setSharpeningLevel(qreal value);

    qreal denoisingLevel() const;
    void setDenoisingLevel(qreal value);

    ColorFilter colorFilter() const;
    void setColorFilter(ColorFilter filter);
    bool isColorFilterSupported(ColorFilter filter) const;

protected:
    ~QCameraImageProcessing();

private:
    friend class QCamera;
    friend class QCameraPrivate;
    QCameraImageProcessing(QCamera *camera);

    Q_DISABLE_COPY(QCameraImageProcessing)
    Q_DECLARE_PRIVATE(QCameraImageProcessing)
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
    QCameraImageProcessingPrivate *d_ptr_deprecated;
#endif
};

QT_END_NAMESPACE

Q_DECLARE_METATYPE(QCameraImageProcessing::WhiteBalanceMode)
Q_DECLARE_METATYPE(QCameraImageProcessing::ColorFilter)

Q_MEDIA_ENUM_DEBUG(QCameraImageProcessing, WhiteBalanceMode)
Q_MEDIA_ENUM_DEBUG(QCameraImageProcessing, ColorFilter)

#endif  // QCAMERAIMAGEPROCESSING_H
