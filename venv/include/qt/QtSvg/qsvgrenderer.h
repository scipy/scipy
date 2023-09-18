/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt SVG module of the Qt Toolkit.
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

#ifndef QSVGRENDERER_H
#define QSVGRENDERER_H

#include <QtGui/qmatrix.h>

#ifndef QT_NO_SVGRENDERER

#include <QtCore/qobject.h>
#include <QtCore/qsize.h>
#include <QtCore/qrect.h>
#include <QtCore/qxmlstream.h>
#include <QtSvg/qtsvgglobal.h>

QT_BEGIN_NAMESPACE


class QSvgRendererPrivate;
class QPainter;
class QByteArray;
class QTransform;

class Q_SVG_EXPORT QSvgRenderer : public QObject
{
    Q_OBJECT

    Q_PROPERTY(QRectF viewBox READ viewBoxF WRITE setViewBox)
    Q_PROPERTY(int framesPerSecond READ framesPerSecond WRITE setFramesPerSecond)
    Q_PROPERTY(int currentFrame READ currentFrame WRITE setCurrentFrame)
    Q_PROPERTY(Qt::AspectRatioMode aspectRatioMode READ aspectRatioMode WRITE setAspectRatioMode)
public:
    QSvgRenderer(QObject *parent = nullptr);
    QSvgRenderer(const QString &filename, QObject *parent = nullptr);
    QSvgRenderer(const QByteArray &contents, QObject *parent = nullptr);
    QSvgRenderer(QXmlStreamReader *contents, QObject *parent = nullptr);
    ~QSvgRenderer();

    bool isValid() const;

    QSize defaultSize() const;

    QRect viewBox() const;
    QRectF viewBoxF() const;
    void setViewBox(const QRect &viewbox);
    void setViewBox(const QRectF &viewbox);

    Qt::AspectRatioMode aspectRatioMode() const;
    void setAspectRatioMode(Qt::AspectRatioMode mode);

    bool animated() const;
    int framesPerSecond() const;
    void setFramesPerSecond(int num);
    int currentFrame() const;
    void setCurrentFrame(int);
    int animationDuration() const;//in seconds

    QRectF boundsOnElement(const QString &id) const;
    bool elementExists(const QString &id) const;
#if QT_DEPRECATED_SINCE(5, 15)
    QT_DEPRECATED_X("Use transformForElement()")
    QMatrix matrixForElement(const QString &id) const;
#endif // QT_DEPRECATED_SINCE(5, 15)
    QTransform transformForElement(const QString &id) const;

public Q_SLOTS:
    bool load(const QString &filename);
    bool load(const QByteArray &contents);
    bool load(QXmlStreamReader *contents);
    void render(QPainter *p);
    void render(QPainter *p, const QRectF &bounds);

    void render(QPainter *p, const QString &elementId,
                const QRectF &bounds=QRectF());

Q_SIGNALS:
    void repaintNeeded();

private:
    Q_DECLARE_PRIVATE(QSvgRenderer)
};

QT_END_NAMESPACE

#endif // QT_NO_SVGRENDERER
#endif // QSVGRENDERER_H
