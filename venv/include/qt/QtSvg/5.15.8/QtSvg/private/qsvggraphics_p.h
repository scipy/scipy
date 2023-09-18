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

#ifndef QSVGGRAPHICS_P_H
#define QSVGGRAPHICS_P_H

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

#include "qsvgnode_p.h"
#include "qtsvgglobal_p.h"

#include "QtGui/qpainterpath.h"
#include "QtGui/qimage.h"
#include "QtGui/qtextlayout.h"
#include "QtGui/qtextoption.h"
#include "QtCore/qstack.h"

QT_BEGIN_NAMESPACE

class QTextCharFormat;

class Q_SVG_PRIVATE_EXPORT QSvgAnimation : public QSvgNode
{
public:
    void draw(QPainter *p, QSvgExtraStates &states) override;
    Type type() const override;
};

class Q_SVG_PRIVATE_EXPORT QSvgArc : public QSvgNode
{
public:
    QSvgArc(QSvgNode *parent, const QPainterPath &path);
    void draw(QPainter *p, QSvgExtraStates &states) override;
    Type type() const override;
    QRectF bounds(QPainter *p, QSvgExtraStates &states) const override;
private:
    QPainterPath m_path;
};

class Q_SVG_PRIVATE_EXPORT QSvgEllipse : public QSvgNode
{
public:
    QSvgEllipse(QSvgNode *parent, const QRectF &rect);
    void draw(QPainter *p, QSvgExtraStates &states) override;
    Type type() const override;
    QRectF bounds(QPainter *p, QSvgExtraStates &states) const override;
private:
    QRectF m_bounds;
};

class Q_SVG_PRIVATE_EXPORT QSvgCircle : public QSvgEllipse
{
public:
    QSvgCircle(QSvgNode *parent, const QRectF &rect) : QSvgEllipse(parent, rect) { }
    Type type() const override;
};

class Q_SVG_PRIVATE_EXPORT QSvgImage : public QSvgNode
{
public:
    QSvgImage(QSvgNode *parent, const QImage &image,
              const QRectF &bounds);
    void draw(QPainter *p, QSvgExtraStates &states) override;
    Type type() const override;
    QRectF bounds(QPainter *p, QSvgExtraStates &states) const override;
private:
    QImage m_image;
    QRectF m_bounds;
};

class Q_SVG_PRIVATE_EXPORT QSvgLine : public QSvgNode
{
public:
    QSvgLine(QSvgNode *parent, const QLineF &line);
    void draw(QPainter *p, QSvgExtraStates &states) override;
    Type type() const override;
    QRectF bounds(QPainter *p, QSvgExtraStates &states) const override;
private:
    QLineF m_line;
};

class Q_SVG_PRIVATE_EXPORT QSvgPath : public QSvgNode
{
public:
    QSvgPath(QSvgNode *parent, const QPainterPath &qpath);
    void draw(QPainter *p, QSvgExtraStates &states) override;
    Type type() const override;
    QRectF bounds(QPainter *p, QSvgExtraStates &states) const override;

    QPainterPath *qpath() {
        return &m_path;
    }
private:
    QPainterPath m_path;
};

class Q_SVG_PRIVATE_EXPORT QSvgPolygon : public QSvgNode
{
public:
    QSvgPolygon(QSvgNode *parent, const QPolygonF &poly);
    void draw(QPainter *p, QSvgExtraStates &states) override;
    Type type() const override;
    QRectF bounds(QPainter *p, QSvgExtraStates &states) const override;
private:
    QPolygonF m_poly;
};

class Q_SVG_PRIVATE_EXPORT QSvgPolyline : public QSvgNode
{
public:
    QSvgPolyline(QSvgNode *parent, const QPolygonF &poly);
    void draw(QPainter *p, QSvgExtraStates &states) override;
    Type type() const override;
    QRectF bounds(QPainter *p, QSvgExtraStates &states) const override;
private:
    QPolygonF m_poly;
};

class Q_SVG_PRIVATE_EXPORT QSvgRect : public QSvgNode
{
public:
    QSvgRect(QSvgNode *paren, const QRectF &rect, int rx=0, int ry=0);
    Type type() const override;
    void draw(QPainter *p, QSvgExtraStates &states) override;
    QRectF bounds(QPainter *p, QSvgExtraStates &states) const override;
private:
    QRectF m_rect;
    int m_rx, m_ry;
};

class  QSvgTspan;

class Q_SVG_PRIVATE_EXPORT QSvgText : public QSvgNode
{
public:
    enum WhitespaceMode
    {
        Default,
        Preserve
    };

    QSvgText(QSvgNode *parent, const QPointF &coord);
    ~QSvgText();
    void setTextArea(const QSizeF &size);

    void draw(QPainter *p, QSvgExtraStates &states) override;
    Type type() const override;

    void addTspan(QSvgTspan *tspan) {m_tspans.append(tspan);}
    void addText(const QString &text);
    void addLineBreak() {m_tspans.append(LINEBREAK);}
    void setWhitespaceMode(WhitespaceMode mode) {m_mode = mode;}

    //QRectF bounds(QPainter *p, QSvgExtraStates &states) const override;
private:
    static QSvgTspan * const LINEBREAK;

    QPointF m_coord;

    // 'm_tspans' is also used to store characters outside tspans and line breaks.
    // If a 'm_tspan' item is null, it indicates a line break.
    QVector<QSvgTspan *> m_tspans;

    Type m_type;
    QSizeF m_size;
    WhitespaceMode m_mode;
};

class Q_SVG_PRIVATE_EXPORT QSvgTspan : public QSvgNode
{
public:
    // tspans are also used to store normal text, so the 'isProperTspan' is used to separate text from tspan.
    QSvgTspan(QSvgNode *parent, bool isProperTspan = true)
        : QSvgNode(parent), m_mode(QSvgText::Default), m_isTspan(isProperTspan)
    {
    }
    ~QSvgTspan() { };
    Type type() const override { return TSPAN; }
    void draw(QPainter *, QSvgExtraStates &) override { Q_ASSERT(!"Tspans should be drawn through QSvgText::draw()."); }
    void addText(const QString &text) {m_text += text;}
    const QString &text() const {return m_text;}
    bool isTspan() const {return m_isTspan;}
    void setWhitespaceMode(QSvgText::WhitespaceMode mode) {m_mode = mode;}
    QSvgText::WhitespaceMode whitespaceMode() const {return m_mode;}
private:
    QString m_text;
    QSvgText::WhitespaceMode m_mode;
    bool m_isTspan;
};

class QSvgUse : public QSvgNode
{
public:
    QSvgUse(const QPointF &start, QSvgNode *parent, QSvgNode *link);
    QSvgUse(const QPointF &start, QSvgNode *parent, const QString &linkId)
        : QSvgUse(start, parent, nullptr)
    { m_linkId = linkId; }
    void draw(QPainter *p, QSvgExtraStates &states) override;
    Type type() const override;
    QRectF bounds(QPainter *p, QSvgExtraStates &states) const override;
    bool isResolved() const { return m_link != nullptr; }
    QString linkId() const { return m_linkId; }
    void setLink(QSvgNode *link) { m_link = link; }

private:
    QSvgNode *m_link;
    QPointF   m_start;
    QString   m_linkId;
    mutable bool m_recursing;
};

class QSvgVideo : public QSvgNode
{
public:
    void draw(QPainter *p, QSvgExtraStates &states) override;
    Type type() const override;
};

QT_END_NAMESPACE

#endif // QSVGGRAPHICS_P_H
