/****************************************************************************
**
** Copyright (C) 2017 The Qt Company Ltd.
** Contact: http://www.qt.io/licensing/
**
** This file is part of the Qt Quick Templates 2 module of the Qt Toolkit.
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

#ifndef QQUICKPALETTE_P_H
#define QQUICKPALETTE_P_H

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

#include <QtGui/qcolor.h>
#include <QtGui/qpalette.h>
#include <QtQuickTemplates2/private/qtquicktemplates2global_p.h>

QT_BEGIN_NAMESPACE

class Q_QUICKTEMPLATES2_PRIVATE_EXPORT QQuickPalette
{
    Q_GADGET
    Q_PROPERTY(QColor alternateBase READ alternateBase WRITE setAlternateBase RESET resetAlternateBase FINAL)
    Q_PROPERTY(QColor base READ base WRITE setBase RESET resetBase FINAL)
    Q_PROPERTY(QColor brightText READ brightText WRITE setBrightText RESET resetBrightText FINAL)
    Q_PROPERTY(QColor button READ button WRITE setButton RESET resetButton FINAL)
    Q_PROPERTY(QColor buttonText READ buttonText WRITE setButtonText RESET resetButtonText FINAL)
    Q_PROPERTY(QColor dark READ dark WRITE setDark RESET resetDark FINAL)
    Q_PROPERTY(QColor highlight READ highlight WRITE setHighlight RESET resetHighlight FINAL)
    Q_PROPERTY(QColor highlightedText READ highlightedText WRITE setHighlightedText RESET resetHighlightedText FINAL)
    Q_PROPERTY(QColor light READ light WRITE setLight RESET resetLight FINAL)
    Q_PROPERTY(QColor link READ link WRITE setLink RESET resetLink FINAL)
    Q_PROPERTY(QColor linkVisited READ linkVisited WRITE setLinkVisited RESET resetLinkVisited FINAL)
    Q_PROPERTY(QColor mid READ mid WRITE setMid RESET resetMid FINAL)
    Q_PROPERTY(QColor midlight READ midlight WRITE setMidlight RESET resetMidlight FINAL)
    Q_PROPERTY(QColor shadow READ shadow WRITE setShadow RESET resetShadow FINAL)
    Q_PROPERTY(QColor text READ text WRITE setText RESET resetText FINAL)
    Q_PROPERTY(QColor toolTipBase READ toolTipBase WRITE setToolTipBase RESET resetToolTipBase FINAL)
    Q_PROPERTY(QColor toolTipText READ toolTipText WRITE setToolTipText RESET resetToolTipText FINAL)
    Q_PROPERTY(QColor window READ window WRITE setWindow RESET resetWindow FINAL)
    Q_PROPERTY(QColor windowText READ windowText WRITE setWindowText RESET resetWindowText FINAL)

public:
    QColor alternateBase() const;
    void setAlternateBase(const QColor &color);
    void resetAlternateBase();

    QColor base() const;
    void setBase(const QColor &color);
    void resetBase();

    QColor brightText() const;
    void setBrightText(const QColor &color);
    void resetBrightText();

    QColor button() const;
    void setButton(const QColor &color);
    void resetButton();

    QColor buttonText() const;
    void setButtonText(const QColor &color);
    void resetButtonText();

    QColor dark() const;
    void setDark(const QColor &color);
    void resetDark();

    QColor highlight() const;
    void setHighlight(const QColor &color);
    void resetHighlight();

    QColor highlightedText() const;
    void setHighlightedText(const QColor &color);
    void resetHighlightedText();

    QColor light() const;
    void setLight(const QColor &color);
    void resetLight();

    QColor link() const;
    void setLink(const QColor &color);
    void resetLink();

    QColor linkVisited() const;
    void setLinkVisited(const QColor &color);
    void resetLinkVisited();

    QColor mid() const;
    void setMid(const QColor &color);
    void resetMid();

    QColor midlight() const;
    void setMidlight(const QColor &color);
    void resetMidlight();

    QColor shadow() const;
    void setShadow(const QColor &color);
    void resetShadow();

    QColor text() const;
    void setText(const QColor &color);
    void resetText();

    QColor toolTipBase() const;
    void setToolTipBase(const QColor &color);
    void resetToolTipBase();

    QColor toolTipText() const;
    void setToolTipText(const QColor &color);
    void resetToolTipText();

    QColor window() const;
    void setWindow(const QColor &color);
    void resetWindow();

    QColor windowText() const;
    void setWindowText(const QColor &color);
    void resetWindowText();

private:
    QPalette v;
};

QT_END_NAMESPACE

#endif // QQUICKPALETTE_P_H
