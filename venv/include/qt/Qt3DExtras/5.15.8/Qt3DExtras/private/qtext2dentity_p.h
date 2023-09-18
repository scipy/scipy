/****************************************************************************
**
** Copyright (C) 2016 Klaralvdalens Datakonsult AB (KDAB).
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

#ifndef QT3DEXTRAS_QTEXT2DENTITY_P_H
#define QT3DEXTRAS_QTEXT2DENTITY_P_H

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API.  It exists for the convenience
// of other Qt classes.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.
//

#include <Qt3DCore/private/qentity_p.h>
#include <Qt3DExtras/private/distancefieldtextrenderer_p.h>
#include <Qt3DExtras/private/qdistancefieldglyphcache_p.h>
#include <QFont>

QT_BEGIN_NAMESPACE

namespace Qt3DCore {
class QScene;
}

namespace Qt3DRender {
class QGeometryRenderer;
class QGeometry;
class QMaterial;
class QAttribute;
class QBuffer;
}

namespace Qt3DExtras {

class QText2DMaterial;
class QText2DEntity;

class QText2DEntityPrivate : public Qt3DCore::QEntityPrivate
{
public:
    QText2DEntityPrivate();
    ~QText2DEntityPrivate();

    Q_DECLARE_PUBLIC(QText2DEntity)

    // keep track of the glyphs currently being displayed,
    // to guarantee proper glyph ref-counting in the
    // QDistanceFieldGlyphCache
    QVector<QGlyphRun> m_currentGlyphRuns;
    QDistanceFieldGlyphCache *m_glyphCache;

    void setScene(Qt3DCore::QScene *scene) override;

    QFont m_font;
    QFont m_scaledFont; // ignore point or pixel size, set to default value

    QColor m_color;
    QString m_text;
    float m_width;
    float m_height;

    QVector<DistanceFieldTextRenderer*> m_renderers;

    float computeActualScale() const;

    void setCurrentGlyphRuns(const QVector<QGlyphRun> &runs);
    void clearCurrentGlyphRuns();
    void update();

    struct CacheEntry
    {
        QDistanceFieldGlyphCache *glyphCache = nullptr;
        int count = 0;
    };

    static QHash<Qt3DCore::QScene *, CacheEntry> m_glyphCacheInstances;
};

} // namespace Qt3DExtras

QT_END_NAMESPACE

#endif // QT3DEXTRAS_QTEXT2DENTITY_P_H
