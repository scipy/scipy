/****************************************************************************
**
** Copyright (C) 2014 Klaralvdalens Datakonsult AB (KDAB).
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt3D module of the Qt Toolkit.
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

#ifndef QT3DRENDER_QSCENE_IMPORTER_P_H
#define QT3DRENDER_QSCENE_IMPORTER_P_H

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

#include <QObject>
#include <QStringList>
#include <QLoggingCategory>
#include <QUrl>
#include <private/qt3drender_global_p.h>

QT_BEGIN_NAMESPACE

namespace Qt3DCore {
class QEntity;
}

namespace Qt3DRender {

Q_DECLARE_LOGGING_CATEGORY(SceneParsers)

class Q_3DRENDERSHARED_PRIVATE_EXPORT QSceneImporter : public QObject
{
    Q_OBJECT
    Q_PROPERTY(ParserStatus status READ status NOTIFY statusChanged)
    Q_PROPERTY(QStringList errors READ errors NOTIFY errorsChanged)

public:
    enum ParserStatus {
        Empty,
        Loading,
        Loaded,
        Error
    };
    Q_ENUM(ParserStatus) // LCOV_EXCL_LINE

    QSceneImporter();
    virtual ~QSceneImporter();

    virtual void setSource(const QUrl &source) = 0;
    virtual void setData(const QByteArray& data, const QString &basePath) = 0;
    virtual bool areFileTypesSupported(const QStringList &extensions) const = 0;
    virtual Qt3DCore::QEntity *scene(const QString &id = QString()) = 0;
    virtual Qt3DCore::QEntity *node(const QString &id) = 0;

    ParserStatus status() const;
    QStringList errors() const;

Q_SIGNALS:
    void statusChanged(ParserStatus status);
    void errorsChanged(const QStringList &errors);

protected:
    void setStatus(ParserStatus status);
    void logError(const QString &error);
    void logInfo(const QString &info);

private:
    ParserStatus m_status;
    QStringList m_errors;
};

} // namespace Qt3DRender

QT_END_NAMESPACE

#endif // QT3DRENDER_QSCENE_IMPORTER_P_H
