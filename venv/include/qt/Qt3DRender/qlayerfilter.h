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

#ifndef QT3DRENDER_QLAYERFILTER_H
#define QT3DRENDER_QLAYERFILTER_H

#include <Qt3DRender/qframegraphnode.h>
#include <QtCore/QStringList>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {

class QLayer;
class QLayerFilterPrivate;

class Q_3DRENDERSHARED_EXPORT QLayerFilter : public QFrameGraphNode
{
    Q_OBJECT
    Q_PROPERTY(FilterMode filterMode READ filterMode WRITE setFilterMode NOTIFY filterModeChanged)
public:
    enum FilterMode
    {
        AcceptAnyMatchingLayers = 0,
        AcceptAllMatchingLayers,
        DiscardAnyMatchingLayers,
        DiscardAllMatchingLayers,
    };
    Q_ENUM(FilterMode) // LOVC_EXLC_LINE

    explicit QLayerFilter(Qt3DCore::QNode *parent = nullptr);
    ~QLayerFilter();

    void addLayer(QLayer *layer);
    void removeLayer(QLayer *layer);
    QVector<QLayer *> layers() const;

    FilterMode filterMode() const;
    void setFilterMode(FilterMode filterMode);

Q_SIGNALS:
    void filterModeChanged(FilterMode filterMode);


protected:
    explicit QLayerFilter(QLayerFilterPrivate &dd, Qt3DCore::QNode *parent = nullptr);

private:
    Q_DECLARE_PRIVATE(QLayerFilter)
    Qt3DCore::QNodeCreatedChangeBasePtr createNodeCreationChange() const override;
};

} // namespace Qt3DRender

QT_END_NAMESPACE

#endif // QT3DRENDER_QLAYERFILTER_H
