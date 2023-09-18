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

#ifndef QSVGSTRUCTURE_P_H
#define QSVGSTRUCTURE_P_H

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

#include "QtCore/qlist.h"
#include "QtCore/qhash.h"

QT_BEGIN_NAMESPACE

class QSvgTinyDocument;
class QSvgNode;
class QPainter;
class QSvgDefs;

class Q_SVG_PRIVATE_EXPORT QSvgStructureNode : public QSvgNode
{
public:
    QSvgStructureNode(QSvgNode *parent);
    ~QSvgStructureNode();
    QSvgNode *scopeNode(const QString &id) const;
    void addChild(QSvgNode *child, const QString &id);
    QRectF bounds(QPainter *p, QSvgExtraStates &states) const override;
    QSvgNode *previousSiblingNode(QSvgNode *n) const;
    QList<QSvgNode*> renderers() const { return m_renderers; }
protected:
    QList<QSvgNode*>          m_renderers;
    QHash<QString, QSvgNode*> m_scope;
    QList<QSvgStructureNode*> m_linkedScopes;
    mutable bool              m_recursing = false;
};

class Q_SVG_PRIVATE_EXPORT QSvgG : public QSvgStructureNode
{
public:
    QSvgG(QSvgNode *parent);
    void draw(QPainter *p, QSvgExtraStates &states) override;
    Type type() const override;
};

class Q_SVG_PRIVATE_EXPORT QSvgDefs : public QSvgStructureNode
{
public:
    QSvgDefs(QSvgNode *parent);
    void draw(QPainter *p, QSvgExtraStates &states) override;
    Type type() const override;
};

class Q_SVG_PRIVATE_EXPORT QSvgSwitch : public QSvgStructureNode
{
public:
    QSvgSwitch(QSvgNode *parent);
    void draw(QPainter *p, QSvgExtraStates &states) override;
    Type type() const override;
private:
    void init();
private:
    QString m_systemLanguage;
    QString m_systemLanguagePrefix;
};

QT_END_NAMESPACE

#endif // QSVGSTRUCTURE_P_H
