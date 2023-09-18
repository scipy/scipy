/****************************************************************************
**
** Copyright (C) 2015 The Qt Company Ltd.
** Contact: http://www.qt.io/licensing/
**
** This file is part of the QtWebView module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:LGPL3$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see http://www.qt.io/terms-conditions. For further
** information use the contact form at http://www.qt.io/contact-us.
**
** GNU Lesser General Public License Usage
** Alternatively, this file may be used under the terms of the GNU Lesser
** General Public License version 3 as published by the Free Software
** Foundation and appearing in the file LICENSE.LGPLv3 included in the
** packaging of this file. Please review the following information to
** ensure the GNU Lesser General Public License version 3 requirements
** will be met: https://www.gnu.org/licenses/lgpl.html.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 2.0 or later as published by the Free
** Software Foundation and appearing in the file LICENSE.GPL included in
** the packaging of this file. Please review the following information to
** ensure the GNU General Public License version 2.0 requirements will be
** met: http://www.gnu.org/licenses/gpl-2.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef QWEBVIEWINTERFACE_H
#define QWEBVIEWINTERFACE_H

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

#include <QtWebView/qwebview_global.h>

#include <QtCore/qobject.h>
#include <QtCore/qurl.h>
#include <QtGui/qimage.h>


QT_BEGIN_NAMESPACE

class QJSValue;

class QWebViewInterface
{
public:
    virtual ~QWebViewInterface() {}
    virtual QString httpUserAgent() const = 0;
    virtual void setHttpUserAgent(const QString &httpUserAgent) = 0;
    virtual QUrl url() const = 0;
    virtual void setUrl(const QUrl &url) = 0;
    virtual bool canGoBack() const = 0;
    virtual bool canGoForward() const = 0;
    virtual QString title() const = 0;
    virtual int loadProgress() const = 0;
    virtual bool isLoading() const = 0;

    // Q_SLOTS
    virtual void goBack() = 0;
    virtual void goForward() = 0;
    virtual void stop() = 0;
    virtual void reload() = 0;
    virtual void loadHtml(const QString &html, const QUrl &baseUrl = QUrl()) = 0;

    virtual void runJavaScriptPrivate(const QString &script,
                                      int callbackId) = 0;
};

QT_END_NAMESPACE

#endif // QWEBVIEWINTERFACE_H
