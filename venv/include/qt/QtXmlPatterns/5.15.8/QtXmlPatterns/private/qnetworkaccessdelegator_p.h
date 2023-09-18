/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtXmlPatterns module of the Qt Toolkit.
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

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API.  It exists purely as an
// implementation detail.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.

#ifndef QNETWORKACCESSDELEGATOR_P_H
#define QNETWORKACCESSDELEGATOR_P_H

#include <QObject>
#include <QPointer>
#include <QSharedData>

QT_BEGIN_NAMESPACE

class QNetworkAccessManager;
class QUrl;

namespace QPatternist
{
    /**
     * @short A value based class that hands out QNetworkAccessManager
     * appropriately for variable bindings.
     *
     * NetworkAccessDelegator is an indirection mechanism which takes care of
     * the fact that QIODevice instances are injected into the data model by
     * having them represented using a URI. managerFor() returns appropriately
     * the right network access manager depending on whether the URI refers to
     * a variable, or to something else.
     *
     * The constructor take a pointer to two NetworkAccessDelegator instances.
     * First is a generic one, the second is the one which handles variable
     * bindings. managerFor() then returns the appropriate one depending on the
     * URI.
     *
     * @since 4.5
     * @see AccelTreeResourceLoader::load()
     * @author Frans Englich <frans.englich@nokia.com>
     */
    class Q_AUTOTEST_EXPORT NetworkAccessDelegator : public QObject
                                                   , public QSharedData
    {
    public:
        typedef QExplicitlySharedDataPointer<NetworkAccessDelegator> Ptr;
        NetworkAccessDelegator(QNetworkAccessManager *const genericManager,
                               QNetworkAccessManager *const variableURIManager);

        QNetworkAccessManager *managerFor(const QUrl &uri);

        QPointer<QNetworkAccessManager> m_genericManager;
        QPointer<QNetworkAccessManager> m_variableURIManager;
    private:
        Q_DISABLE_COPY(NetworkAccessDelegator)
    };
}

QT_END_NAMESPACE
#endif
