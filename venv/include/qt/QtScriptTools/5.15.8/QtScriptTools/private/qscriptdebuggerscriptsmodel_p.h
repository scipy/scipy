/****************************************************************************
**
** Copyright (C) 2018 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtSCriptTools module of the Qt Toolkit.
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

#ifndef QSCRIPTDEBUGGERSCRIPTSMODEL_P_H
#define QSCRIPTDEBUGGERSCRIPTSMODEL_P_H

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

#include <QtCore/qabstractitemmodel.h>
#include <QtCore/qpair.h>

#include "qscriptscriptdata_p.h"

QT_BEGIN_NAMESPACE

class QScriptDebuggerScriptsModelPrivate;
class Q_AUTOTEST_EXPORT QScriptDebuggerScriptsModel
    : public QAbstractItemModel
{
    Q_OBJECT
public:
    QScriptDebuggerScriptsModel(QObject *parent = 0);
    ~QScriptDebuggerScriptsModel();

    void removeScript(qint64 id);
    void addScript(qint64 id, const QScriptScriptData &data);
    void addExtraScriptInfo(
        qint64 id, const QMap<QString, int> &functionsInfo,
        const QSet<int> &executableLineNumbers);
    void commit();

    QScriptScriptData scriptData(qint64 id) const;
    QScriptScriptMap scripts() const;
    qint64 resolveScript(const QString &fileName) const;
    QSet<int> executableLineNumbers(qint64 scriptId) const;

    QModelIndex indexFromScriptId(qint64 id) const;
    qint64 scriptIdFromIndex(const QModelIndex &index) const;
    QPair<QString, int> scriptFunctionInfoFromIndex(const QModelIndex &index) const;

    QModelIndex index(int row, int column, const QModelIndex &parent = QModelIndex()) const;
    QModelIndex parent(const QModelIndex &child) const;
    int columnCount(const QModelIndex &parent = QModelIndex()) const;
    int rowCount(const QModelIndex &parent = QModelIndex()) const;
    QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const;

private:
    Q_DECLARE_PRIVATE(QScriptDebuggerScriptsModel)
    Q_DISABLE_COPY(QScriptDebuggerScriptsModel)
};

QT_END_NAMESPACE

#endif
