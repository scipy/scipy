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

#ifndef QWEBVIEW_P_H
#define QWEBVIEW_P_H

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

#include "qabstractwebview_p.h"
#include "qwebviewinterface_p.h"
#include "qnativeviewcontroller_p.h"
#include <QtCore/qobject.h>
#include <QtCore/qurl.h>
#include <QtGui/qimage.h>
#include <QtQml/qjsvalue.h>

QT_BEGIN_NAMESPACE

class QWebViewLoadRequestPrivate;

class Q_WEBVIEW_EXPORT QWebView
        : public QObject
        , public QWebViewInterface
        , public QNativeViewController
{
    Q_OBJECT
public:
    enum LoadStatus { // Changes here needs to be done in QQuickWebView as well
        LoadStartedStatus,
        LoadStoppedStatus,
        LoadSucceededStatus,
        LoadFailedStatus
    };

    explicit QWebView(QObject *p = 0);
    ~QWebView() Q_DECL_OVERRIDE;

    QString httpUserAgent() const Q_DECL_OVERRIDE;
    void setHttpUserAgent(const QString &httpUserAgent) Q_DECL_OVERRIDE;
    QUrl url() const Q_DECL_OVERRIDE;
    void setUrl(const QUrl &url) Q_DECL_OVERRIDE;
    bool canGoBack() const Q_DECL_OVERRIDE;
    bool canGoForward() const Q_DECL_OVERRIDE;
    QString title() const Q_DECL_OVERRIDE;
    int loadProgress() const Q_DECL_OVERRIDE;
    bool isLoading() const Q_DECL_OVERRIDE;

    void setParentView(QObject *view) Q_DECL_OVERRIDE;
    QObject *parentView() const Q_DECL_OVERRIDE;
    void setGeometry(const QRect &geometry) Q_DECL_OVERRIDE;
    void setVisibility(QWindow::Visibility visibility) Q_DECL_OVERRIDE;
    void setVisible(bool visible) Q_DECL_OVERRIDE;
    void setFocus(bool focus) Q_DECL_OVERRIDE;

public Q_SLOTS:
    void goBack() Q_DECL_OVERRIDE;
    void goForward() Q_DECL_OVERRIDE;
    void reload() Q_DECL_OVERRIDE;
    void stop() Q_DECL_OVERRIDE;
    void loadHtml(const QString &html, const QUrl &baseUrl = QUrl()) Q_DECL_OVERRIDE;

Q_SIGNALS:
    void titleChanged();
    void urlChanged();
    void loadingChanged(const QWebViewLoadRequestPrivate &loadRequest);
    void loadProgressChanged();
    void javaScriptResult(int id, const QVariant &result);
    void requestFocus(bool focus);
    void httpUserAgentChanged();

protected:
    void init() Q_DECL_OVERRIDE;
    void runJavaScriptPrivate(const QString &script,
                              int callbackId) Q_DECL_OVERRIDE;

private Q_SLOTS:
    void onTitleChanged(const QString &title);
    void onUrlChanged(const QUrl &url);
    void onLoadProgressChanged(int progress);
    void onLoadingChanged(const QWebViewLoadRequestPrivate &loadRequest);
    void onHttpUserAgentChanged(const QString &httpUserAgent);

private:
    friend class QQuickViewController;
    friend class QQuickWebView;

    QAbstractWebView *d;

    // provisional data
    int m_progress;
    QString m_title;
    QUrl m_url;
    mutable QString m_httpUserAgent;
};

QT_END_NAMESPACE

#endif // QWEBVIEW_P_H
