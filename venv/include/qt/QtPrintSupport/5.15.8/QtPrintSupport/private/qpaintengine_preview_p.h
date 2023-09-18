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

#ifndef QPAINTENGINE_PREVIEW_P_H
#define QPAINTENGINE_PREVIEW_P_H

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API.  It exists for the convenience
// of QPreviewPrinter and friends.  This header file may change from
// version to version without notice, or even be removed.
//
// We mean it.
//
//

#include <QtPrintSupport/private/qtprintsupportglobal_p.h>
#include <QtGui/qpaintengine.h>
#include <QtPrintSupport/qprintengine.h>

QT_REQUIRE_CONFIG(printpreviewwidget);

QT_BEGIN_NAMESPACE

class QPreviewPaintEnginePrivate;

class QPreviewPaintEngine : public QPaintEngine, public QPrintEngine
{
    Q_DECLARE_PRIVATE(QPreviewPaintEngine)
public:
    QPreviewPaintEngine();
    ~QPreviewPaintEngine();

    bool begin(QPaintDevice *dev) override;
    bool end() override;

    void updateState(const QPaintEngineState &state) override;

    void drawPath(const QPainterPath &path) override;
    void drawPolygon(const QPointF *points, int pointCount, PolygonDrawMode mode) override;
    void drawTextItem(const QPointF &p, const QTextItem &textItem) override;

    void drawPixmap(const QRectF &r, const QPixmap &pm, const QRectF &sr) override;
    void drawTiledPixmap(const QRectF &r, const QPixmap &pm, const QPointF &p) override;

    QList<const QPicture *> pages();

    QPaintEngine::Type type() const override { return Picture; }

    void setProxyEngines(QPrintEngine *printEngine, QPaintEngine *paintEngine);

    void setProperty(PrintEnginePropertyKey key, const QVariant &value) override;
    QVariant property(PrintEnginePropertyKey key) const override;

    bool newPage() override;
    bool abort() override;

    int metric(QPaintDevice::PaintDeviceMetric) const override;

    QPrinter::PrinterState printerState() const override;
};

QT_END_NAMESPACE

#endif
