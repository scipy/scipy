/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt Data Visualization module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:GPL$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 3 or (at your option) any later version
** approved by the KDE Free Qt Foundation. The licenses are as published by
** the Free Software Foundation and appearing in the file LICENSE.GPL3
** included in the packaging of this file. Please review the following
** information to ensure the GNU General Public License requirements will
** be met: https://www.gnu.org/licenses/gpl-3.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

//
//  W A R N I N G
//  -------------
//
// This file is not part of the QtDataVisualization API.  It exists purely as an
// implementation detail.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.

#ifndef QCUSTOM3DITEM_P_H
#define QCUSTOM3DITEM_P_H

#include "datavisualizationglobal_p.h"
#include "qcustom3ditem.h"

QT_BEGIN_NAMESPACE_DATAVISUALIZATION

struct QCustomItemDirtyBitField {
    bool textureDirty               : 1;
    bool meshDirty                  : 1;
    bool positionDirty              : 1;
    bool scalingDirty               : 1;
    bool rotationDirty              : 1;
    bool visibleDirty               : 1;
    bool shadowCastingDirty         : 1;

    QCustomItemDirtyBitField()
        : textureDirty(false),
          meshDirty(false),
          positionDirty(false),
          scalingDirty(false),
          rotationDirty(false),
          visibleDirty(false),
          shadowCastingDirty(false)
    {
    }
};

class QCustom3DItemPrivate : public QObject
{
    Q_OBJECT
public:
    QCustom3DItemPrivate(QCustom3DItem *q);
    QCustom3DItemPrivate(QCustom3DItem *q, const QString &meshFile, const QVector3D &position,
                         const QVector3D &scaling, const QQuaternion &rotation);
    virtual ~QCustom3DItemPrivate();

    QImage textureImage();
    void clearTextureImage();
    void resetDirtyBits();

public:
    QCustom3DItem *q_ptr;
    QImage m_textureImage;
    QString m_textureFile;
    QString m_meshFile;
    QVector3D m_position;
    bool m_positionAbsolute;
    QVector3D m_scaling;
    bool m_scalingAbsolute;
    QQuaternion m_rotation;
    bool m_visible;
    bool m_shadowCasting;

    bool m_isLabelItem;
    bool m_isVolumeItem;

    QCustomItemDirtyBitField m_dirtyBits;

Q_SIGNALS:
    void needUpdate();

private:
    QCustom3DItemPrivate(QCustom3DItemPrivate *d);

    friend class QCustom3DItem;
};

QT_END_NAMESPACE_DATAVISUALIZATION

#endif
