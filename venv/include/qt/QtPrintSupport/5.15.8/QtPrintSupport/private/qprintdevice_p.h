/****************************************************************************
**
** Copyright (C) 2014 John Layt <jlayt@kde.org>
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtPrintSupport module of the Qt Toolkit.
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

#ifndef QPRINTDEVICE_H
#define QPRINTDEVICE_H

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API.  It exists for the convenience
// of internal files.  This header file may change from version to version
// without notice, or even be removed.
//
// We mean it.
//

#include <QtPrintSupport/private/qtprintsupportglobal_p.h>
#include "private/qprint_p.h"

#include <QtCore/qsharedpointer.h>
#include <QtGui/qpagelayout.h>

QT_BEGIN_NAMESPACE

#ifndef QT_NO_PRINTER

class QPlatformPrintDevice;
class QMarginsF;
class QMimeType;
class QDebug;

class Q_PRINTSUPPORT_EXPORT QPrintDevice
{
public:

    QPrintDevice();
    QPrintDevice(const QString & id);
    QPrintDevice(const QPrintDevice &other);
    ~QPrintDevice();

    QPrintDevice &operator=(const QPrintDevice &other);
    QPrintDevice &operator=(QPrintDevice &&other) { swap(other); return *this; }

    void swap(QPrintDevice &other) { d.swap(other.d); }

    bool operator==(const QPrintDevice &other) const;

    QString id() const;
    QString name() const;
    QString location() const;
    QString makeAndModel() const;

    bool isValid() const;
    bool isDefault() const;
    bool isRemote() const;

    QPrint::DeviceState state() const;

    bool isValidPageLayout(const QPageLayout &layout, int resolution) const;

    bool supportsMultipleCopies() const;
    bool supportsCollateCopies() const;

    QPageSize defaultPageSize() const;
    QList<QPageSize> supportedPageSizes() const;

    QPageSize supportedPageSize(const QPageSize &pageSize) const;
    QPageSize supportedPageSize(QPageSize::PageSizeId pageSizeId) const;
    QPageSize supportedPageSize(const QString &pageName) const;
    QPageSize supportedPageSize(const QSize &pointSize) const;
    QPageSize supportedPageSize(const QSizeF &size, QPageSize::Unit units = QPageSize::Point) const;

    bool supportsCustomPageSizes() const;

    QSize minimumPhysicalPageSize() const;
    QSize maximumPhysicalPageSize() const;

    QMarginsF printableMargins(const QPageSize &pageSize, QPageLayout::Orientation orientation, int resolution) const;

    int defaultResolution() const;
    QList<int> supportedResolutions() const;

    QPrint::InputSlot defaultInputSlot() const;
    QVector<QPrint::InputSlot> supportedInputSlots() const;

    QPrint::OutputBin defaultOutputBin() const;
    QVector<QPrint::OutputBin> supportedOutputBins() const;

    QPrint::DuplexMode defaultDuplexMode() const;
    QVector<QPrint::DuplexMode> supportedDuplexModes() const;

    QPrint::ColorMode defaultColorMode() const;
    QVector<QPrint::ColorMode> supportedColorModes() const;

    enum PrintDevicePropertyKey {
        PDPK_CustomBase = 0xff00
    };

    QVariant property(PrintDevicePropertyKey key) const;
    bool setProperty(PrintDevicePropertyKey key, const QVariant &value);
    bool isFeatureAvailable(PrintDevicePropertyKey key, const QVariant &params) const;

#if QT_CONFIG(mimetype)
    QList<QMimeType> supportedMimeTypes() const;
#endif

#  ifndef QT_NO_DEBUG_STREAM
    void format(QDebug debug) const;
#  endif

private:
    friend class QPlatformPrinterSupport;
    friend class QPlatformPrintDevice;
    QPrintDevice(QPlatformPrintDevice *dd);
    QSharedPointer<QPlatformPrintDevice> d;
};

Q_DECLARE_SHARED(QPrintDevice)

#  ifndef QT_NO_DEBUG_STREAM
Q_PRINTSUPPORT_EXPORT QDebug operator<<(QDebug debug, const QPrintDevice &);
#  endif
#endif // QT_NO_PRINTER

QT_END_NAMESPACE

#endif // QPLATFORMPRINTDEVICE_H
