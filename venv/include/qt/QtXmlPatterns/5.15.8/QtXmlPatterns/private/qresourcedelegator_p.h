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

#ifndef QPatternist_ResourceDelegator_p_H
#define QPatternist_ResourceDelegator_p_H

#include <QSet>
#include <QUrl>

#include <private/qdeviceresourceloader_p.h>

QT_BEGIN_NAMESPACE

namespace QPatternist
{
    /**
     * @short Delegates to another ResourceLoader, but in case a URI is in an
     * exception list, it delegates to a different loader.
     *
     * This is used for handling device variables, since when a device variable
     * is rebound, a resource loader needs to carry that binding, while the
     * resource loader for the other query remains as is.
     *
     * @since 4.5
     * @author Frans Englich <frans.englich@nokia.com>
     */
    class ResourceDelegator : public DeviceResourceLoader
    {
    public:
        ResourceDelegator(const QSet<QUrl> &needsOverride,
                          const ResourceLoader::Ptr &parentLoader,
                          const ResourceLoader::Ptr &forDeviceLoader) : m_needsOverride(needsOverride)
                                                                      , m_parentLoader(parentLoader)
                                                                      , m_forDeviceLoader(forDeviceLoader)

        {
            Q_ASSERT(m_parentLoader);
        }

        bool isUnparsedTextAvailable(const QUrl &uri,
                                     const QString &encoding) override;
        ItemType::Ptr announceUnparsedText(const QUrl &uri) override;
        Item openUnparsedText(const QUrl &uri,
                              const QString &encoding,
                              const ReportContext::Ptr &context,
                              const SourceLocationReflection *const where) override;
        Item openDocument(const QUrl &uri,
                          const ReportContext::Ptr &context) override;
        SequenceType::Ptr announceDocument(const QUrl &uri, const Usage usageHint) override;
        bool isDocumentAvailable(const QUrl &uri) override;
        Item::Iterator::Ptr openCollection(const QUrl &uri) override;
        SequenceType::Ptr announceCollection(const QUrl &uri) override;

        /**
         * Returns the union of the deviceURIs() that ResourceDelegator's two
         * resource loaders has.
         */
        QSet<QUrl> deviceURIs() const override;

    private:
        const QSet<QUrl> m_needsOverride;
        const ResourceLoader::Ptr m_parentLoader;
        const ResourceDelegator::Ptr m_forDeviceLoader;
    };
}

QT_END_NAMESPACE

#endif
