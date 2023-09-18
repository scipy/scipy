/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtGui module of the Qt Toolkit.
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

#ifndef QCUPS_P_H
#define QCUPS_P_H

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

#include <QtPrintSupport/private/qtprintsupportglobal_p.h>
#include <QtPrintSupport/private/qprint_p.h>
#include "QtCore/qstring.h"
#include "QtCore/qstringlist.h"
#include "QtPrintSupport/qprinter.h"
#include "QtCore/qdatetime.h"

QT_REQUIRE_CONFIG(cups);

QT_BEGIN_NAMESPACE

class QPrintDevice;

// HACK! Define these here temporarily so they can be used in the dialogs
// without a circular reference to QCupsPrintEngine in the plugin.
// Move back to qcupsprintengine_p.h in the plugin once all usage
// removed from the dialogs.
#define PPK_CupsOptions QPrintEngine::PrintEnginePropertyKey(0xfe00)

#define PDPK_PpdFile          QPrintDevice::PrintDevicePropertyKey(QPrintDevice::PDPK_CustomBase)
#define PDPK_PpdOption        QPrintDevice::PrintDevicePropertyKey(QPrintDevice::PDPK_CustomBase + 1)
#define PDPK_CupsJobPriority  QPrintDevice::PrintDevicePropertyKey(QPrintDevice::PDPK_CustomBase + 2)
#define PDPK_CupsJobSheets    QPrintDevice::PrintDevicePropertyKey(QPrintDevice::PDPK_CustomBase + 3)
#define PDPK_CupsJobBilling   QPrintDevice::PrintDevicePropertyKey(QPrintDevice::PDPK_CustomBase + 4)
#define PDPK_CupsJobHoldUntil QPrintDevice::PrintDevicePropertyKey(QPrintDevice::PDPK_CustomBase + 5)
#define PDPK_PpdChoiceIsInstallableConflict QPrintDevice::PrintDevicePropertyKey(QPrintDevice::PDPK_CustomBase + 6)

class Q_PRINTSUPPORT_EXPORT QCUPSSupport
{
public:
    // Enum for values of job-hold-until option
    enum JobHoldUntil {
        NoHold = 0,  //CUPS Default
        Indefinite,
        DayTime,
        Night,
        SecondShift,
        ThirdShift,
        Weekend,
        SpecificTime
    };

    // Enum for valid banner pages
    enum BannerPage {
        NoBanner = 0,  //CUPS Default 'none'
        Standard,
        Unclassified,
        Confidential,
        Classified,
        Secret,
        TopSecret
    };

    // Enum for valid page set
    enum PageSet {
        AllPages = 0,  //CUPS Default
        OddPages,
        EvenPages
    };

    // Enum for valid number of pages per sheet
    enum PagesPerSheet {
        OnePagePerSheet = 0,
        TwoPagesPerSheet,
        FourPagesPerSheet,
        SixPagesPerSheet,
        NinePagesPerSheet,
        SixteenPagesPerSheet
    };

    // Enum for valid layouts of pages per sheet
    enum PagesPerSheetLayout {
        LeftToRightTopToBottom = 0,
        LeftToRightBottomToTop,
        RightToLeftTopToBottom,
        RightToLeftBottomToTop,
        BottomToTopLeftToRight,
        BottomToTopRightToLeft,
        TopToBottomLeftToRight,
        TopToBottomRightToLeft
    };

    static void setCupsOption(QPrinter *printer, const QString &option, const QString &value);
    static void clearCupsOption(QPrinter *printer, const QString &option);
    static void clearCupsOptions(QPrinter *printer);

    static void setJobHold(QPrinter *printer, const JobHoldUntil jobHold = NoHold, const QTime &holdUntilTime = QTime());
    static void setJobBilling(QPrinter *printer, const QString &jobBilling = QString());
    static void setJobPriority(QPrinter *printer, int priority = 50);
    static void setBannerPages(QPrinter *printer, const BannerPage startBannerPage, const BannerPage endBannerPage);
    static void setPageSet(QPrinter *printer, const PageSet pageSet);
    static void setPagesPerSheetLayout(QPrinter *printer, const PagesPerSheet pagesPerSheet,
                                       const PagesPerSheetLayout pagesPerSheetLayout);
    static void setPageRange(QPrinter *printer, int pageFrom, int pageTo);
    static void setPageRange(QPrinter *printer, const QString &pageRange);

    struct JobSheets
    {
        JobSheets(BannerPage s = NoBanner, BannerPage e = NoBanner)
         : startBannerPage(s), endBannerPage(e) {}

        BannerPage startBannerPage;
        BannerPage endBannerPage;
    };
    static JobSheets parseJobSheets(const QString &jobSheets);

    struct JobHoldUntilWithTime
    {
        JobHoldUntilWithTime(JobHoldUntil jh = NoHold, const QTime &t = QTime())
            : jobHold(jh), time(t) {}

        JobHoldUntil jobHold;
        QTime time;
    };
    static JobHoldUntilWithTime parseJobHoldUntil(const QString &jobHoldUntil);

    static ppd_option_t *findPpdOption(const char *optionName, QPrintDevice *printDevice);
};
Q_DECLARE_TYPEINFO(QCUPSSupport::JobHoldUntil,        Q_PRIMITIVE_TYPE);
Q_DECLARE_TYPEINFO(QCUPSSupport::BannerPage,          Q_PRIMITIVE_TYPE);
Q_DECLARE_TYPEINFO(QCUPSSupport::PageSet,             Q_PRIMITIVE_TYPE);
Q_DECLARE_TYPEINFO(QCUPSSupport::PagesPerSheetLayout, Q_PRIMITIVE_TYPE);
Q_DECLARE_TYPEINFO(QCUPSSupport::PagesPerSheet,       Q_PRIMITIVE_TYPE);

QT_END_NAMESPACE

Q_DECLARE_METATYPE(QCUPSSupport::JobHoldUntil)
Q_DECLARE_METATYPE(QCUPSSupport::BannerPage)
Q_DECLARE_METATYPE(QCUPSSupport::PageSet)
Q_DECLARE_METATYPE(QCUPSSupport::PagesPerSheetLayout)
Q_DECLARE_METATYPE(QCUPSSupport::PagesPerSheet)

#endif
