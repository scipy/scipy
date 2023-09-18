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

#ifndef Patternist_AccelTreeResourceLoader_H
#define Patternist_AccelTreeResourceLoader_H

#include <QtCore/QEventLoop>
#include <QtNetwork/QNetworkReply>

#include <QAbstractXmlReceiver>
#include <private/qacceltree_p.h>
#include <private/qacceltreebuilder_p.h>
#include <private/qdeviceresourceloader_p.h>
#include <private/qnamepool_p.h>
#include <private/qnetworkaccessdelegator_p.h>
#include <private/qreportcontext_p.h>

QT_BEGIN_NAMESPACE

class QIODevice;

namespace QPatternist
{
    /**
     * @short Handles requests for documents, and instantiates
     * them as AccelTree instances.
     *
     * @author Frans Englich <frans.englich@nokia.com>
     */
    class Q_AUTOTEST_EXPORT AccelTreeResourceLoader : public DeviceResourceLoader
    {
    public:
        /**
         * Describes the behaviour of the resource loader in case of an
         * error.
         */
        enum ErrorHandling
        {
            FailOnError,        ///< The resource loader will report the error via the report context.
            ContinueOnError     ///< The resource loader will report no error and return an empty QNetworkReply.
        };

        /**
         * AccelTreeResourceLoader does not own @p context.
         */
        AccelTreeResourceLoader(const NamePool::Ptr &np,
                                const NetworkAccessDelegator::Ptr &networkDelegator, AccelTreeBuilder<true>::Features = AccelTreeBuilder<true>::NoneFeature);

        virtual Item openDocument(const QUrl &uri,
                                  const ReportContext::Ptr &context);
        virtual Item openDocument(QIODevice *source, const QUrl &documentUri,
                                  const ReportContext::Ptr &context);
        virtual SequenceType::Ptr announceDocument(const QUrl &uri, const Usage usageHint);
        virtual bool isDocumentAvailable(const QUrl &uri);

        virtual bool isUnparsedTextAvailable(const QUrl &uri,
                                             const QString &encoding);

        virtual Item openUnparsedText(const QUrl &uri,
                                      const QString &encoding,
                                      const ReportContext::Ptr &context,
                                      const SourceLocationReflection *const where);

        /**
         * @short Helper function that do NetworkAccessDelegator::get(), but
         * does it blocked.
         *
         * The returned QNetworkReply has emitted QNetworkReply::finished().
         *
         * The caller owns the return QIODevice instance.
         *
         * @p context may be @c null or valid. If @c null, no error reporting
         * is done and @c null is returned.
         *
         * @see NetworkAccessDelegator
         */
        static QNetworkReply *load(const QUrl &uri,
                                   QNetworkAccessManager *const networkManager,
                                   const ReportContext::Ptr &context, ErrorHandling handling = FailOnError);

        /**
         * @overload
         */
        static QNetworkReply *load(const QUrl &uri,
                                   const NetworkAccessDelegator::Ptr &networkDelegator,
                                   const ReportContext::Ptr &context, ErrorHandling handling = FailOnError);

        /**
         * @short Returns the URIs this AccelTreeResourceLoader has loaded
         * which are for devices through variable bindings.
         */
        virtual QSet<QUrl> deviceURIs() const;

        virtual void clear(const QUrl &uri);

    private:
        static bool streamToReceiver(QIODevice *const dev,
                                     AccelTreeBuilder<true> *const receiver,
                                     const NamePool::Ptr &np,
                                     const ReportContext::Ptr &context,
                                     const QUrl &uri);
        bool retrieveDocument(const QUrl &uri,
                              const ReportContext::Ptr &context);
        bool retrieveDocument(QIODevice *source, const QUrl &documentUri,
                              const ReportContext::Ptr &context);
        /**
         * If @p context is @c null, no error reporting should be done.
         */
        bool retrieveUnparsedText(const QUrl &uri,
                                  const QString &encoding,
                                  const ReportContext::Ptr &context,
                                  const SourceLocationReflection *const where);

        QHash<QUrl, AccelTree::Ptr>             m_loadedDocuments;
        const NamePool::Ptr                     m_namePool;
        const NetworkAccessDelegator::Ptr       m_networkAccessDelegator;
        QHash<QPair<QUrl, QString>, QString>    m_unparsedTexts;
        AccelTreeBuilder<true>::Features        m_features;
    };
}

QT_END_NAMESPACE

#endif
