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

#ifndef QPRINTER_P_H
#define QPRINTER_P_H

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

#ifndef QT_NO_PRINTER

#include "QtPrintSupport/qprinter.h"
#include "QtPrintSupport/qprinterinfo.h"
#include "QtPrintSupport/qprintengine.h"
#include "QtCore/qpointer.h"
#include "QtCore/qset.h"

#include <limits.h>

QT_BEGIN_NAMESPACE

class QPrintEngine;
class QPreviewPaintEngine;
class QPicture;

class Q_PRINTSUPPORT_EXPORT QPrinterPrivate
{
    Q_DECLARE_PUBLIC(QPrinter)
public:
    QPrinterPrivate(QPrinter *printer)
        : pdfVersion(QPrinter::PdfVersion_1_4),
          printEngine(nullptr),
          paintEngine(nullptr),
          realPrintEngine(nullptr),
          realPaintEngine(nullptr),
#if QT_CONFIG(printpreviewwidget)
          previewEngine(nullptr),
#endif
          q_ptr(printer),
          printRange(QPrinter::AllPages),
          use_default_engine(true),
          validPrinter(false)
    {
    }

    ~QPrinterPrivate() {

    }

    static QPrinterPrivate *get(QPrinter *printer) {
        return printer->d_ptr.get();
    }

    void init(const QPrinterInfo &printer, QPrinter::PrinterMode mode);

    QPrinterInfo findValidPrinter(const QPrinterInfo &printer = QPrinterInfo());
    void initEngines(QPrinter::OutputFormat format, const QPrinterInfo &printer);
    void changeEngines(QPrinter::OutputFormat format, const QPrinterInfo &printer);
#if QT_CONFIG(printpreviewwidget)
    QList<const QPicture *> previewPages() const;
    void setPreviewMode(bool);
#endif

    void setProperty(QPrintEngine::PrintEnginePropertyKey key, const QVariant &value);

    QPrinter::PrinterMode printerMode;
    QPrinter::OutputFormat outputFormat;
    QPrinter::PdfVersion pdfVersion;
    QPrintEngine *printEngine;
    QPaintEngine *paintEngine;

    QPrintEngine *realPrintEngine;
    QPaintEngine *realPaintEngine;
#if QT_CONFIG(printpreviewwidget)
    QPreviewPaintEngine *previewEngine;
#endif

    QPrinter *q_ptr;

    QPrinter::PrintRange printRange;

    uint use_default_engine : 1;
    uint had_default_engines : 1;

    uint validPrinter : 1;
    uint hasCustomPageMargins : 1;

    // Used to remember which properties have been manually set by the user.
    QSet<QPrintEngine::PrintEnginePropertyKey> m_properties;
};

QT_END_NAMESPACE

#endif // QT_NO_PRINTER

#endif // QPRINTER_P_H
