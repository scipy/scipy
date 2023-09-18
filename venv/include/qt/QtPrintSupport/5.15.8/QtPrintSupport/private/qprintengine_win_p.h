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

#ifndef QPRINTENGINE_WIN_P_H
#define QPRINTENGINE_WIN_P_H

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

#include <QtPrintSupport/private/qtprintsupportglobal_p.h>

#ifndef QT_NO_PRINTER

#include <QtGui/qpaintengine.h>
#include <QtGui/qpagelayout.h>
#include <QtPrintSupport/QPrintEngine>
#include <QtPrintSupport/QPrinter>
#include <private/qpaintengine_alpha_p.h>
#include <private/qprintdevice_p.h>
#include <QtCore/qt_windows.h>

QT_BEGIN_NAMESPACE

class QWin32PrintEnginePrivate;
class QPrinterPrivate;
class QPainterState;

class Q_PRINTSUPPORT_EXPORT QWin32PrintEngine : public QAlphaPaintEngine, public QPrintEngine
{
    Q_DECLARE_PRIVATE(QWin32PrintEngine)
public:
    QWin32PrintEngine(QPrinter::PrinterMode mode, const QString &deviceId);

    // override QWin32PaintEngine
    bool begin(QPaintDevice *dev);
    bool end();

    void updateState(const QPaintEngineState &state);

    void updateMatrix(const QTransform &matrix);
    void updateClipPath(const QPainterPath &clip, Qt::ClipOperation op);

    void drawPath(const QPainterPath &path);
    void drawPolygon(const QPointF *points, int pointCount, PolygonDrawMode mode);
    void drawTextItem(const QPointF &p, const QTextItem &textItem);

    void drawPixmap(const QRectF &r, const QPixmap &pm, const QRectF &sr);
    void drawTiledPixmap(const QRectF &r, const QPixmap &pm, const QPointF &p);
    void setProperty(PrintEnginePropertyKey key, const QVariant &value);
    QVariant property(PrintEnginePropertyKey key) const;

    bool newPage();
    bool abort();
    int metric(QPaintDevice::PaintDeviceMetric) const;

    QPrinter::PrinterState printerState() const;

    QPaintEngine::Type type() const { return Windows; }

    HDC getDC() const;
    void releaseDC(HDC) const;

    /* Used by print/page setup dialogs */
    void setGlobalDevMode(HGLOBAL globalDevNames, HGLOBAL globalDevMode);
    HGLOBAL *createGlobalDevNames();
    HGLOBAL globalDevMode();

private:
    friend class QPrintDialog;
    friend class QPageSetupDialog;
};

class QWin32PrintEnginePrivate : public QAlphaPaintEnginePrivate
{
    Q_DECLARE_PUBLIC(QWin32PrintEngine)
public:
    QWin32PrintEnginePrivate() :
        hPrinter(0),
        globalDevMode(0),
        devMode(0),
        pInfo(0),
        hMem(0),
        hdc(0),
        ownsDevMode(false),
        mode(QPrinter::ScreenResolution),
        state(QPrinter::Idle),
        resolution(0),
        m_pageLayout(QPageLayout(QPageSize(QPageSize::A4), QPageLayout::Portrait, QMarginsF(0, 0, 0, 0))),
        stretch_x(1), stretch_y(1), origin_x(0), origin_y(0),
        dpi_x(96), dpi_y(96), dpi_display(96),
        num_copies(1),
        printToFile(false),
        reinit(false),
        complex_xform(false), has_pen(false), has_brush(false), has_custom_paper_size(false),
        embed_fonts(true),
        txop(0 /* QTransform::TxNone */)
    {
    }

    ~QWin32PrintEnginePrivate();


    /* Initializes the printer data based on the current printer name. This
       function creates a DEVMODE struct, HDC and a printer handle. If these
       structures are already in use, they are freed using release
    */
    void initialize();

    /* Initializes data in the print engine whenever the HDC has been renewed
    */
    void initHDC();

    /* Releases all the handles the printer currently holds, HDC, DEVMODE,
       etc and resets the corresponding members to 0. */
    void release();

    /* Resets the DC with changes in devmode. If the printer is active
       this function only sets the reinit variable to true so it
       is handled in the next begin or newpage. */
    void doReinit();

    bool resetDC();

    void strokePath(const QPainterPath &path, const QColor &color);
    void fillPath(const QPainterPath &path, const QColor &color);

    void composeGdiPath(const QPainterPath &path);
    void fillPath_dev(const QPainterPath &path, const QColor &color);
    void strokePath_dev(const QPainterPath &path, const QColor &color, qreal width);

    void setPageSize(const QPageSize &pageSize);
    void updatePageLayout();
    void updateMetrics();
    void debugMetrics() const;

    // Windows GDI printer references.
    HANDLE hPrinter;

    HGLOBAL globalDevMode;
    DEVMODE *devMode;
    PRINTER_INFO_2 *pInfo;
    HGLOBAL hMem;

    HDC hdc;

    // True if devMode was allocated separately from pInfo.
    bool ownsDevMode;

    QPrinter::PrinterMode mode;

    // Print Device
    QPrintDevice m_printDevice;

    // Document info
    QString docName;
    QString m_creator;
    QString fileName;

    QPrinter::PrinterState state;
    int resolution;

    // Page Layout
    QPageLayout m_pageLayout;

    // Page metrics cache
    QRect m_paintRectPixels;
    QSize m_paintSizeMM;

    // Windows painting
    qreal stretch_x;
    qreal stretch_y;
    int origin_x;
    int origin_y;

    int dpi_x;
    int dpi_y;
    int dpi_display;
    int num_copies;

    uint printToFile : 1;
    uint reinit : 1;

    uint complex_xform : 1;
    uint has_pen : 1;
    uint has_brush : 1;
    uint has_custom_paper_size : 1;
    uint embed_fonts : 1;

    uint txop;

    QColor brush_color;
    QPen pen;
    QColor pen_color;
    QSizeF paper_size;  // In points

    QTransform painterMatrix;
    QTransform matrix;
};

QT_END_NAMESPACE

#endif // QT_NO_PRINTER

#endif // QPRINTENGINE_WIN_P_H
