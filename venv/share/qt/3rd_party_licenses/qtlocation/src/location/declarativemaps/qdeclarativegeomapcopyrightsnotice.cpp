/****************************************************************************
**
** Copyright (C) 2014 Aaron McCarthy <mccarthy.aaron@gmail.com>
** Copyright (C) 2015 The Qt Company Ltd.
** Contact: http://www.qt.io/licensing/
**
** This file is part of the QtLocation module of the Qt Toolkit.
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

#include "qdeclarativegeomapcopyrightsnotice_p.h"

#include <QtGui/QTextDocument>
#include <QtGui/QAbstractTextDocumentLayout>
#include <QtGui/QPainter>
#include <QtQuick/private/qquickanchors_p.h>
#include <QtQuick/private/qquickanchors_p_p.h>
#include <QtLocation/private/qdeclarativegeomap_p.h>
#include <QtQuick/private/qquickpainteditem_p.h>

QT_BEGIN_NAMESPACE

class QDeclarativeGeoMapCopyrightNoticePrivate: public QQuickPaintedItemPrivate
{
    Q_DECLARE_PUBLIC(QDeclarativeGeoMapCopyrightNotice)
public:
    virtual void setVisible(bool visible);
};

/*!
    \qmltype MapCopyrightNotice
    \instantiates QDeclarativeGeoMapCopyrightNotice
    \inqmlmodule QtLocation
    \ingroup qml-QtLocation5-maps
    \since QtLocation 5.9

    \brief The MapCopyrightNotice item displays the current valid
    copyright notice for a Map element.

    This object can be used to place an additional copyright notices
    programmatically.

    Note that declaring a MapCopyrightNotice inside a QtLocation::Map element
    is not possible, like for any other QQuickItem.

    The release of this API with Qt 5.9 is a Technology Preview.
*/

/*!
    \qmlproperty Map QtLocation::MapCopyrightNotice::mapSource

    This property holds the current map source providing the copyright data shown
    in this notice.
    In order to let the MapCopyrightNotice display a copyright, this property must
    be set, as it is the only data source for this element.
*/

/*!
    \qmlproperty string QtLocation::MapCopyrightNotice::styleSheet

    This property holds the current css2.1 style sheet used to style the copyright notice, if in HTML form.

    Example:
    \code
    MapCopyrightNotice {
        mapSource: myMap
        styleSheet: "body { color : green; font-family: \"Lucida\"; font-size: 8px} a{ font-size: 8px; color:#A62900}"
    }
    \endcode
*/

QDeclarativeGeoMapCopyrightNotice::QDeclarativeGeoMapCopyrightNotice(QQuickItem *parent)
:   QQuickPaintedItem(parent), m_copyrightsHtml(0), m_copyrightsVisible(true), m_mapSource(0),
    m_userDefinedStyleSheet(false)
{
    // If this item is constructed inside the map, automatically anchor it where it always used to be.
    if (qobject_cast<QDeclarativeGeoMap *>(parent))
        anchorToBottomLeft();
}

QDeclarativeGeoMapCopyrightNotice::~QDeclarativeGeoMapCopyrightNotice()
{
    setMapSource(nullptr);
}

void QDeclarativeGeoMapCopyrightNotice::anchorToBottomLeft()
{
    if (!parent())
        return;
    QQuickAnchors *anchors = property("anchors").value<QQuickAnchors *>();
    if (anchors) {
        anchors->setLeft(QQuickAnchorLine(qobject_cast<QQuickItem *>(parent()), QQuickAnchors::LeftAnchor));
        anchors->setBottom(QQuickAnchorLine(qobject_cast<QQuickItem *>(parent()), QQuickAnchors::BottomAnchor));
    }
}

void QDeclarativeGeoMapCopyrightNotice::setMapSource(QDeclarativeGeoMap *map)
{
    if (m_mapSource == map)
        return;

    if (m_mapSource) {
        // disconnect this object from current map source
        m_mapSource->detachCopyrightNotice(copyrightsVisible());
        m_mapSource->disconnect(this);
        m_mapSource->m_map->disconnect(this);
        if (m_copyrightsHtml)
            m_copyrightsHtml->clear();
        m_copyrightsImage = QImage();
        m_mapSource = nullptr;
    }

    if (map) {
        m_mapSource = map;
        m_mapSource->attachCopyrightNotice(copyrightsVisible());
        connect(this, &QDeclarativeGeoMapCopyrightNotice::copyrightsVisibleChanged,
                mapSource(), &QDeclarativeGeoMap::onAttachedCopyrightNoticeVisibilityChanged);

        // First update the copyright. Only Image will do here, no need to store HTML right away.
        if (m_mapSource->m_copyrights && !m_mapSource->m_copyrights->m_copyrightsImage.isNull())
            m_copyrightsImage = m_mapSource->m_copyrights->m_copyrightsImage;

        connect(mapSource(), SIGNAL(copyrightsChanged(QImage)),
                this, SLOT(copyrightsChanged(QImage)));
        connect(mapSource(), SIGNAL(copyrightsChanged(QString)),
                this, SLOT(copyrightsChanged(QString)));

        if (m_mapSource->m_map)
            connectMap();
        else
            connect(mapSource(), &QDeclarativeGeoMap::mapReadyChanged, this, &QDeclarativeGeoMapCopyrightNotice::connectMap);
    }
}

void QDeclarativeGeoMapCopyrightNotice::connectMap()
{
    connect(m_mapSource->m_map, SIGNAL(copyrightsStyleSheetChanged(QString)),
            this, SLOT(onCopyrightsStyleSheetChanged(QString)));
    connect(this, SIGNAL(linkActivated(QString)),
            mapSource(), SIGNAL(copyrightLinkActivated(QString)));

    onCopyrightsStyleSheetChanged(m_mapSource->m_map->copyrightsStyleSheet());

    update();
    emit mapSourceChanged();
}

QDeclarativeGeoMap *QDeclarativeGeoMapCopyrightNotice::mapSource()
{
    return m_mapSource.data();
}

QString QDeclarativeGeoMapCopyrightNotice::styleSheet() const
{
    return m_styleSheet;
}

void QDeclarativeGeoMapCopyrightNotice::setStyleSheet(const QString &styleSheet)
{
    m_userDefinedStyleSheet = true;

    if (styleSheet == m_styleSheet)
        return;

    m_styleSheet = styleSheet;
    if (!m_html.isEmpty() && m_copyrightsHtml) {
        delete m_copyrightsHtml;
        createCopyright();
#if QT_CONFIG(texthtmlparser)
        m_copyrightsHtml->setHtml(m_html);
#else
        m_copyrightsHtml->setPlainText(m_html);
#endif
    }
    rasterizeHtmlAndUpdate();
    emit styleSheetChanged(m_styleSheet);
}

/*!
    \internal
*/
void QDeclarativeGeoMapCopyrightNotice::paint(QPainter *painter)
{
    painter->drawImage(0, 0, m_copyrightsImage);
}

void QDeclarativeGeoMapCopyrightNotice::mousePressEvent(QMouseEvent *event)
{
    if (m_copyrightsHtml) {
        m_activeAnchor = m_copyrightsHtml->documentLayout()->anchorAt(event->pos());
        if (!m_activeAnchor.isEmpty())
            return;
    }

    QQuickPaintedItem::mousePressEvent(event);
}

void QDeclarativeGeoMapCopyrightNotice::mouseReleaseEvent(QMouseEvent *event)
{
    if (m_copyrightsHtml) {
        QString anchor = m_copyrightsHtml->documentLayout()->anchorAt(event->pos());
        if (anchor == m_activeAnchor && !anchor.isEmpty()) {
            emit linkActivated(anchor);
            m_activeAnchor.clear();
        }
    }
}

void QDeclarativeGeoMapCopyrightNotice::rasterizeHtmlAndUpdate()
{
    if (!m_copyrightsHtml || m_copyrightsHtml->isEmpty())
        return;

    m_copyrightsImage = QImage(m_copyrightsHtml->size().toSize(),
                               QImage::Format_ARGB32_Premultiplied);

    m_copyrightsImage.fill(qPremultiply(QColor(Qt::transparent).rgba()));
    QPainter painter(&m_copyrightsImage);
    QAbstractTextDocumentLayout::PaintContext ctx;
    ctx.palette.setColor(QPalette::Text, QStringLiteral("black"));
    m_copyrightsHtml->documentLayout()->draw(&painter, ctx);

    setImplicitSize(m_copyrightsImage.width(), m_copyrightsImage.height());
    setContentsSize(m_copyrightsImage.size());

    setKeepMouseGrab(true);
    setAcceptedMouseButtons(Qt::LeftButton);

    update();
}

void QDeclarativeGeoMapCopyrightNotice::createCopyright()
{
    m_copyrightsHtml = new QTextDocument(this);
#if QT_CONFIG(cssparser)
    if (!m_styleSheet.isEmpty())
        m_copyrightsHtml->setDefaultStyleSheet(m_styleSheet);
#endif

    // The default 4 makes the copyright too wide and tall.
    m_copyrightsHtml->setDocumentMargin(0);
}

void QDeclarativeGeoMapCopyrightNoticePrivate::setVisible(bool visible)
{
    Q_Q(QDeclarativeGeoMapCopyrightNotice);
    q->m_copyrightsVisible = visible;
    QQuickItemPrivate::setVisible(visible);
}

/*!
    \internal
*/
void QDeclarativeGeoMapCopyrightNotice::setCopyrightsVisible(bool visible)
{
    Q_D(QDeclarativeGeoMapCopyrightNotice);
    if (visible == m_copyrightsVisible)
        return;

    m_copyrightsVisible = visible;
    d->QQuickItemPrivate::setVisible(!m_copyrightsImage.isNull() && visible);
    emit copyrightsVisibleChanged();
}

bool QDeclarativeGeoMapCopyrightNotice::copyrightsVisible() const
{
    return m_copyrightsVisible;
}

/*!
    \internal
*/
void QDeclarativeGeoMapCopyrightNotice::setCopyrightsZ(qreal copyrightsZ)
{
    setZ(copyrightsZ);
    update();
}

/*!
    \internal
*/
void QDeclarativeGeoMapCopyrightNotice::copyrightsChanged(const QImage &copyrightsImage)
{
    Q_D(QDeclarativeGeoMapCopyrightNotice);
    delete m_copyrightsHtml;
    m_copyrightsHtml = 0;

    m_copyrightsImage = copyrightsImage;

    setImplicitSize(m_copyrightsImage.width(), m_copyrightsImage.height());

    setKeepMouseGrab(false);
    setAcceptedMouseButtons(Qt::NoButton);
    d->QQuickItemPrivate::setVisible(m_copyrightsVisible && !m_copyrightsImage.isNull());

    update();
}

void QDeclarativeGeoMapCopyrightNotice::copyrightsChanged(const QString &copyrightsHtml)
{
    Q_D(QDeclarativeGeoMapCopyrightNotice);
    if (copyrightsHtml.isEmpty()) {
        d->QQuickItemPrivate::setVisible(false);
        return;
    } else  {
        d->QQuickItemPrivate::setVisible(m_copyrightsVisible);
    }

    // Divfy, so we can style the background. The extra <span> is a
    // workaround to QTBUG-58838 and should be removed when it gets fixed.
#if QT_CONFIG(texthtmlparser)
    m_html = QStringLiteral("<div id='copyright-root'><span>") + copyrightsHtml + QStringLiteral("</span></div>");
#else
    m_html = copyrightsHtml;
#endif

    if (!m_copyrightsHtml)
        createCopyright();

#if QT_CONFIG(texthtmlparser)
    m_copyrightsHtml->setHtml(m_html);
#else
    m_copyrightsHtml->setPlainText(m_html);
#endif
    rasterizeHtmlAndUpdate();
}

void QDeclarativeGeoMapCopyrightNotice::onCopyrightsStyleSheetChanged(const QString &styleSheet)
{
    if (m_userDefinedStyleSheet || styleSheet == m_styleSheet)
        return;

    m_styleSheet = styleSheet;
    if (!m_html.isEmpty() && m_copyrightsHtml) {
        delete m_copyrightsHtml;
        createCopyright();
#if QT_CONFIG(texthtmlparser)
        m_copyrightsHtml->setHtml(m_html);
#else
        m_copyrightsHtml->setPlainText(m_html);
#endif
    }
    rasterizeHtmlAndUpdate();
    emit styleSheetChanged(m_styleSheet);
}

QT_END_NAMESPACE
