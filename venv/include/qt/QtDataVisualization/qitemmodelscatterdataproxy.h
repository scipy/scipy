/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt Data Visualization module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:GPL$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 3 or (at your option) any later version
** approved by the KDE Free Qt Foundation. The licenses are as published by
** the Free Software Foundation and appearing in the file LICENSE.GPL3
** included in the packaging of this file. Please review the following
** information to ensure the GNU General Public License requirements will
** be met: https://www.gnu.org/licenses/gpl-3.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef QITEMMODELSCATTERDATAPROXY_H
#define QITEMMODELSCATTERDATAPROXY_H

#include <QtDataVisualization/qscatterdataproxy.h>
#include <QtCore/QAbstractItemModel>
#include <QtCore/QString>
#include <QtCore/QRegExp>

QT_BEGIN_NAMESPACE_DATAVISUALIZATION

class QItemModelScatterDataProxyPrivate;

class QT_DATAVISUALIZATION_EXPORT QItemModelScatterDataProxy : public QScatterDataProxy
{
    Q_OBJECT
    Q_PROPERTY(QAbstractItemModel* itemModel READ itemModel WRITE setItemModel NOTIFY itemModelChanged)
    Q_PROPERTY(QString xPosRole READ xPosRole WRITE setXPosRole NOTIFY xPosRoleChanged)
    Q_PROPERTY(QString yPosRole READ yPosRole WRITE setYPosRole NOTIFY yPosRoleChanged)
    Q_PROPERTY(QString zPosRole READ zPosRole WRITE setZPosRole NOTIFY zPosRoleChanged)
    Q_PROPERTY(QString rotationRole READ rotationRole WRITE setRotationRole NOTIFY rotationRoleChanged)
    Q_PROPERTY(QRegExp xPosRolePattern READ xPosRolePattern WRITE setXPosRolePattern NOTIFY xPosRolePatternChanged REVISION 1)
    Q_PROPERTY(QRegExp yPosRolePattern READ yPosRolePattern WRITE setYPosRolePattern NOTIFY yPosRolePatternChanged REVISION 1)
    Q_PROPERTY(QRegExp zPosRolePattern READ zPosRolePattern WRITE setZPosRolePattern NOTIFY zPosRolePatternChanged REVISION 1)
    Q_PROPERTY(QRegExp rotationRolePattern READ rotationRolePattern WRITE setRotationRolePattern NOTIFY rotationRolePatternChanged REVISION 1)
    Q_PROPERTY(QString xPosRoleReplace READ xPosRoleReplace WRITE setXPosRoleReplace NOTIFY xPosRoleReplaceChanged REVISION 1)
    Q_PROPERTY(QString yPosRoleReplace READ yPosRoleReplace WRITE setYPosRoleReplace NOTIFY yPosRoleReplaceChanged REVISION 1)
    Q_PROPERTY(QString zPosRoleReplace READ zPosRoleReplace WRITE setZPosRoleReplace NOTIFY zPosRoleReplaceChanged REVISION 1)
    Q_PROPERTY(QString rotationRoleReplace READ rotationRoleReplace WRITE setRotationRoleReplace NOTIFY rotationRoleReplaceChanged REVISION 1)

public:
    explicit QItemModelScatterDataProxy(QObject *parent = nullptr);
    explicit QItemModelScatterDataProxy(QAbstractItemModel *itemModel, QObject *parent = nullptr);
    explicit QItemModelScatterDataProxy(QAbstractItemModel *itemModel,
                                        const QString &xPosRole, const QString &yPosRole,
                                        const QString &zPosRole, QObject *parent = nullptr);
    explicit QItemModelScatterDataProxy(QAbstractItemModel *itemModel,
                                        const QString &xPosRole, const QString &yPosRole,
                                        const QString &zPosRole, const QString &rotationRole,
                                        QObject *parent = nullptr);
    virtual ~QItemModelScatterDataProxy();

    void setItemModel(QAbstractItemModel *itemModel);
    QAbstractItemModel *itemModel() const;

    void setXPosRole(const QString &role);
    QString xPosRole() const;
    void setYPosRole(const QString &role);
    QString yPosRole() const;
    void setZPosRole(const QString &role);
    QString zPosRole() const;
    void setRotationRole(const QString &role);
    QString rotationRole() const;

    void remap(const QString &xPosRole, const QString &yPosRole, const QString &zPosRole,
               const QString &rotationRole);

    void setXPosRolePattern(const QRegExp &pattern);
    QRegExp xPosRolePattern() const;
    void setYPosRolePattern(const QRegExp &pattern);
    QRegExp yPosRolePattern() const;
    void setZPosRolePattern(const QRegExp &pattern);
    QRegExp zPosRolePattern() const;
    void setRotationRolePattern(const QRegExp &pattern);
    QRegExp rotationRolePattern() const;

    void setXPosRoleReplace(const QString &replace);
    QString xPosRoleReplace() const;
    void setYPosRoleReplace(const QString &replace);
    QString yPosRoleReplace() const;
    void setZPosRoleReplace(const QString &replace);
    QString zPosRoleReplace() const;
    void setRotationRoleReplace(const QString &replace);
    QString rotationRoleReplace() const;

Q_SIGNALS:
    void itemModelChanged(const QAbstractItemModel* itemModel);
    void xPosRoleChanged(const QString &role);
    void yPosRoleChanged(const QString &role);
    void zPosRoleChanged(const QString &role);
    void rotationRoleChanged(const QString &role);
    Q_REVISION(1) void xPosRolePatternChanged(const QRegExp &pattern);
    Q_REVISION(1) void yPosRolePatternChanged(const QRegExp &pattern);
    Q_REVISION(1) void zPosRolePatternChanged(const QRegExp &pattern);
    Q_REVISION(1) void rotationRolePatternChanged(const QRegExp &pattern);
    Q_REVISION(1) void rotationRoleReplaceChanged(const QString &replace);
    Q_REVISION(1) void xPosRoleReplaceChanged(const QString &replace);
    Q_REVISION(1) void yPosRoleReplaceChanged(const QString &replace);
    Q_REVISION(1) void zPosRoleReplaceChanged(const QString &replace);

protected:
    QItemModelScatterDataProxyPrivate *dptr();
    const QItemModelScatterDataProxyPrivate *dptrc() const;

private:
    Q_DISABLE_COPY(QItemModelScatterDataProxy)

    friend class ScatterItemModelHandler;
};

QT_END_NAMESPACE_DATAVISUALIZATION

#endif
