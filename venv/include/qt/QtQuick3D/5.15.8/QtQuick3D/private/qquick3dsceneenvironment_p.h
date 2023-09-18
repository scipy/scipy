/****************************************************************************
**
** Copyright (C) 2019 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of Qt Quick 3D.
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

#ifndef QSSGSCENEENVIRONMENT_H
#define QSSGSCENEENVIRONMENT_H

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

#include <QtCore/QObject>
#include <QtCore/QVector>

#include <QtGui/QColor>

#include <QtQuick3D/private/qquick3dnode_p.h>

#include <QtQml/QQmlListProperty>

#include <QtQuick3D/private/qquick3deffect_p.h>

QT_BEGIN_NAMESPACE

class QQuick3DTexture;
class Q_QUICK3D_EXPORT QQuick3DSceneEnvironment : public QQuick3DObject
{
    Q_OBJECT
    Q_PROPERTY(QQuick3DEnvironmentAAModeValues antialiasingMode READ antialiasingMode WRITE setAntialiasingMode NOTIFY antialiasingModeChanged)
    Q_PROPERTY(QQuick3DEnvironmentAAQualityValues antialiasingQuality READ antialiasingQuality WRITE setAntialiasingQuality NOTIFY antialiasingQualityChanged)

    Q_PROPERTY(bool temporalAAEnabled READ temporalAAEnabled WRITE setTemporalAAEnabled NOTIFY temporalAAEnabledChanged)
    Q_PROPERTY(float temporalAAStrength READ temporalAAStrength WRITE setTemporalAAStrength NOTIFY temporalAAStrengthChanged)
    Q_PROPERTY(QQuick3DEnvironmentBackgroundTypes backgroundMode READ backgroundMode WRITE setBackgroundMode NOTIFY backgroundModeChanged)
    Q_PROPERTY(QColor clearColor READ clearColor WRITE setClearColor NOTIFY clearColorChanged)
    Q_PROPERTY(bool depthTestEnabled READ depthTestEnabled WRITE setDepthTestEnabled NOTIFY depthTestEnabledChanged)
    Q_PROPERTY(bool depthPrePassEnabled READ depthPrePassEnabled WRITE setDepthPrePassEnabled NOTIFY depthPrePassEnabledChanged)

    Q_PROPERTY(float aoStrength READ aoStrength WRITE setAoStrength NOTIFY aoStrengthChanged)
    Q_PROPERTY(float aoDistance READ aoDistance WRITE setAoDistance NOTIFY aoDistanceChanged)
    Q_PROPERTY(float aoSoftness READ aoSoftness WRITE setAoSoftness NOTIFY aoSoftnessChanged)
    Q_PROPERTY(bool aoDither READ aoDither WRITE setAoDither NOTIFY aoDitherChanged)
    Q_PROPERTY(int aoSampleRate READ aoSampleRate WRITE setAoSampleRate NOTIFY aoSampleRateChanged)
    Q_PROPERTY(float aoBias READ aoBias WRITE setAoBias NOTIFY aoBiasChanged)

    Q_PROPERTY(QQuick3DTexture *lightProbe READ lightProbe WRITE setLightProbe NOTIFY lightProbeChanged)
    Q_PROPERTY(float probeBrightness READ probeBrightness WRITE setProbeBrightness NOTIFY probeBrightnessChanged)
    Q_PROPERTY(bool fastImageBasedLightingEnabled READ fastImageBasedLightingEnabled WRITE setFastImageBasedLightingEnabled NOTIFY fastImageBasedLightingEnabledChanged)
    Q_PROPERTY(float probeHorizon READ probeHorizon WRITE setProbeHorizon NOTIFY probeHorizonChanged)
    Q_PROPERTY(float probeFieldOfView READ probeFieldOfView WRITE setProbeFieldOfView NOTIFY probeFieldOfViewChanged)
    Q_PROPERTY(QQmlListProperty<QQuick3DEffect> effects READ effects REVISION 1)

public:

    enum QQuick3DEnvironmentAAModeValues {
        NoAA = 0,
        SSAA,
        MSAA,
        ProgressiveAA
    };
    Q_ENUM(QQuick3DEnvironmentAAModeValues)

    enum QQuick3DEnvironmentAAQualityValues {
        Medium = 2,
        High = 4,
        VeryHigh = 8
    };
    Q_ENUM(QQuick3DEnvironmentAAQualityValues)

    enum QQuick3DEnvironmentBackgroundTypes {
        Transparent = 0,
        Unspecified,
        Color,
        SkyBox
    };
    Q_ENUM(QQuick3DEnvironmentBackgroundTypes)

    explicit QQuick3DSceneEnvironment(QQuick3DObject *parent = nullptr);
    ~QQuick3DSceneEnvironment() override;

    QQuick3DEnvironmentAAModeValues antialiasingMode() const;
    QQuick3DEnvironmentAAQualityValues antialiasingQuality() const;
    bool temporalAAEnabled() const;
    float temporalAAStrength() const;

    QQuick3DEnvironmentBackgroundTypes backgroundMode() const;
    QColor clearColor() const;

    float aoStrength() const;
    float aoDistance() const;
    float aoSoftness() const;
    bool aoDither() const;
    int aoSampleRate() const;
    float aoBias() const;

    QQuick3DTexture *lightProbe() const;
    float probeBrightness() const;
    bool fastImageBasedLightingEnabled() const;
    float probeHorizon() const;
    float probeFieldOfView() const;

    bool depthTestEnabled() const;
    bool depthPrePassEnabled() const;

    QQmlListProperty<QQuick3DEffect> effects();

public Q_SLOTS:
    void setAntialiasingMode(QQuick3DEnvironmentAAModeValues antialiasingMode);
    void setAntialiasingQuality(QQuick3DEnvironmentAAQualityValues antialiasingQuality);
    void setTemporalAAEnabled(bool temporalAAEnabled);
    void setTemporalAAStrength(float strength);

    void setBackgroundMode(QQuick3DEnvironmentBackgroundTypes backgroundMode);
    void setClearColor(const QColor &clearColor);

    void setAoStrength(float aoStrength);
    void setAoDistance(float aoDistance);
    void setAoSoftness(float aoSoftness);
    void setAoDither(bool aoDither);
    void setAoSampleRate(int aoSampleRate);
    void setAoBias(float aoBias);

    void setLightProbe(QQuick3DTexture *lightProbe);
    void setProbeBrightness(float probeBrightness);
    void setFastImageBasedLightingEnabled(bool fastImageBasedLightingEnabled);
    void setProbeHorizon(float probeHorizon);
    void setProbeFieldOfView(float probeFieldOfView);

    void setDepthTestEnabled(bool depthTestEnabled);
    void setDepthPrePassEnabled(bool depthPrePassEnabled);

Q_SIGNALS:
    void antialiasingModeChanged();
    void antialiasingQualityChanged();
    void temporalAAEnabledChanged();
    void temporalAAStrengthChanged();

    void backgroundModeChanged();
    void clearColorChanged();

    void aoStrengthChanged();
    void aoDistanceChanged();
    void aoSoftnessChanged();
    void aoDitherChanged();
    void aoSampleRateChanged();
    void aoBiasChanged();

    void lightProbeChanged();
    void probeBrightnessChanged();
    void fastImageBasedLightingEnabledChanged();
    void probeHorizonChanged();
    void probeFieldOfViewChanged();

    void depthTestEnabledChanged();
    void depthPrePassEnabledChanged();

protected:
    QSSGRenderGraphObject *updateSpatialNode(QSSGRenderGraphObject *node) override;
    void itemChange(ItemChange, const ItemChangeData &) override;

private:
    friend class QQuick3DSceneRenderer;

    QVector<QQuick3DEffect *> m_effects;

    static void qmlAppendEffect(QQmlListProperty<QQuick3DEffect> *list, QQuick3DEffect *effect);
    static QQuick3DEffect *qmlEffectAt(QQmlListProperty<QQuick3DEffect> *list, int index);
    static int qmlEffectsCount(QQmlListProperty<QQuick3DEffect> *list);
    static void qmlClearEffects(QQmlListProperty<QQuick3DEffect> *list);

    void updateSceneManager(const QSharedPointer<QQuick3DSceneManager> &manager);

    QQuick3DEnvironmentAAModeValues m_antialiasingMode = NoAA;
    QQuick3DEnvironmentAAQualityValues m_antialiasingQuality = High;
    bool m_temporalAAEnabled = false;
    float m_temporalAAStrength = 0.3f;

    QQuick3DEnvironmentBackgroundTypes m_backgroundMode = Transparent;
    QColor m_clearColor = Qt::black;

    float m_aoStrength = 0.0f;
    float m_aoDistance = 5.0f;
    float m_aoSoftness = 50.0f;
    bool m_aoDither = false;
    int m_aoSampleRate = 2;
    float m_aoBias = 0.0f;
    QQuick3DTexture *m_lightProbe = nullptr;
    float m_probeBrightness = 100.0f;
    bool m_fastImageBasedLightingEnabled = false;
    float m_probeHorizon = -1.0f;
    float m_probeFieldOfView = 180.0f;

    ConnectionMap m_connections;
    bool m_depthTestEnabled = true;
    bool m_depthPrePassEnabled = false;
};

QT_END_NAMESPACE

#endif // QSSGSCENEENVIRONMENT_H
