/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtScxml module of the Qt Toolkit.
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

#ifndef QSCXMLSTATEMACHINE_H
#define QSCXMLSTATEMACHINE_H

#include <QtScxml/qscxmldatamodel.h>
#include <QtScxml/qscxmlerror.h>
#include <QtScxml/qscxmlevent.h>
#include <QtScxml/qscxmlcompiler.h>
#include <QtScxml/qscxmlinvokableservice.h>

#include <QtCore/qstring.h>
#include <QtCore/qvector.h>
#include <QtCore/qurl.h>
#include <QtCore/qvariant.h>
#include <QtCore/qpointer.h>

#include <functional>

QT_BEGIN_NAMESPACE
class QIODevice;
class QXmlStreamWriter;
class QTextStream;

class QScxmlStateMachinePrivate;
class Q_SCXML_EXPORT QScxmlStateMachine: public QObject
{
    Q_DECLARE_PRIVATE(QScxmlStateMachine)
    Q_OBJECT
    Q_PROPERTY(bool running READ isRunning WRITE setRunning NOTIFY runningChanged)
    Q_PROPERTY(bool initialized READ isInitialized NOTIFY initializedChanged)
    Q_PROPERTY(QScxmlDataModel *dataModel READ dataModel WRITE setDataModel NOTIFY dataModelChanged)
    Q_PROPERTY(QVariantMap initialValues READ initialValues WRITE setInitialValues NOTIFY initialValuesChanged)
    Q_PROPERTY(QVector<QScxmlInvokableService*> invokedServices READ invokedServices NOTIFY invokedServicesChanged)
    Q_PROPERTY(QString sessionId READ sessionId CONSTANT)
    Q_PROPERTY(QString name READ name CONSTANT)
    Q_PROPERTY(bool invoked READ isInvoked CONSTANT)
    Q_PROPERTY(QVector<QScxmlError> parseErrors READ parseErrors CONSTANT)
    Q_PROPERTY(QScxmlCompiler::Loader *loader READ loader WRITE setLoader NOTIFY loaderChanged)
    Q_PROPERTY(QScxmlTableData *tableData READ tableData WRITE setTableData NOTIFY tableDataChanged)

protected:
    explicit QScxmlStateMachine(const QMetaObject *metaObject, QObject *parent = nullptr);
    QScxmlStateMachine(QScxmlStateMachinePrivate &dd, QObject *parent = nullptr);

public:
    static QScxmlStateMachine *fromFile(const QString &fileName);
    static QScxmlStateMachine *fromData(QIODevice *data, const QString &fileName = QString());
    QVector<QScxmlError> parseErrors() const;

    QString sessionId() const;

    bool isInvoked() const;
    bool isInitialized() const;

    void setDataModel(QScxmlDataModel *model);
    QScxmlDataModel *dataModel() const;

    void setLoader(QScxmlCompiler::Loader *loader);
    QScxmlCompiler::Loader *loader() const;

    bool isRunning() const;
    void setRunning(bool running);

    QVariantMap initialValues();
    void setInitialValues(const QVariantMap &initialValues);

    QString name() const;
    Q_INVOKABLE QStringList stateNames(bool compress = true) const;
    Q_INVOKABLE QStringList activeStateNames(bool compress = true) const;
    Q_INVOKABLE bool isActive(const QString &scxmlStateName) const;

    QMetaObject::Connection connectToState(const QString &scxmlStateName,
                                           const QObject *receiver, const char *method,
                                           Qt::ConnectionType type = Qt::AutoConnection);

    // connect state to a QObject slot
    template <typename PointerToMemberFunction>
    inline QMetaObject::Connection connectToState(
            const QString &scxmlStateName,
            const typename QtPrivate::FunctionPointer<PointerToMemberFunction>::Object *receiver,
            PointerToMemberFunction method,
            Qt::ConnectionType type = Qt::AutoConnection)
    {
        typedef QtPrivate::FunctionPointer<PointerToMemberFunction> SlotType;
        return connectToStateImpl(
                    scxmlStateName, receiver, nullptr,
                    new QtPrivate::QSlotObject<PointerToMemberFunction,
                    typename SlotType::Arguments, void>(method),
                    type);
    }

    // connect state to a functor or function pointer (without context)
    template <typename Functor>
    inline typename QtPrivate::QEnableIf<
            !QtPrivate::FunctionPointer<Functor>::IsPointerToMemberFunction &&
            !std::is_same<const char*, Functor>::value, QMetaObject::Connection>::Type
    connectToState(const QString &scxmlStateName, Functor functor,
                   Qt::ConnectionType type = Qt::AutoConnection)
    {
        // Use this as context
        return connectToState(scxmlStateName, this, functor, type);
    }

    // connectToState to a functor or function pointer (with context)
    template <typename Functor>
    inline typename QtPrivate::QEnableIf<
            !QtPrivate::FunctionPointer<Functor>::IsPointerToMemberFunction &&
            !std::is_same<const char*, Functor>::value, QMetaObject::Connection>::Type
    connectToState(const QString &scxmlStateName, const QObject *context, Functor functor,
                   Qt::ConnectionType type = Qt::AutoConnection)
    {
        QtPrivate::QSlotObjectBase *slotObj = new QtPrivate::QFunctorSlotObject<Functor, 1,
                QtPrivate::List<bool>, void>(functor);
        return connectToStateImpl(scxmlStateName, context, reinterpret_cast<void **>(&functor),
                                  slotObj, type);
    }

    static std::function<void(bool)> onEntry(const QObject *receiver, const char *method)
    {
        const QPointer<QObject> receiverPointer(const_cast<QObject *>(receiver));
        return [receiverPointer, method](bool isEnteringState) {
            if (isEnteringState && !receiverPointer.isNull())
                QMetaObject::invokeMethod(const_cast<QObject *>(receiverPointer.data()), method);
        };
    }

    static std::function<void(bool)> onExit(const QObject *receiver, const char *method)
    {
        const QPointer<QObject> receiverPointer(const_cast<QObject *>(receiver));
        return [receiverPointer, method](bool isEnteringState) {
            if (!isEnteringState && !receiverPointer.isNull())
                QMetaObject::invokeMethod(receiverPointer.data(), method);
        };
    }

    template<typename Functor>
    static std::function<void(bool)> onEntry(Functor functor)
    {
        return [functor](bool isEnteringState) {
            if (isEnteringState)
                functor();
        };
    }

    template<typename Functor>
    static std::function<void(bool)> onExit(Functor functor)
    {
        return [functor](bool isEnteringState) {
            if (!isEnteringState)
                functor();
        };
    }

    template<typename PointerToMemberFunction>
    static std::function<void(bool)> onEntry(
            const typename QtPrivate::FunctionPointer<PointerToMemberFunction>::Object *receiver,
            PointerToMemberFunction method)
    {
        typedef typename QtPrivate::FunctionPointer<PointerToMemberFunction>::Object Object;
        const QPointer<Object> receiverPointer(const_cast<Object *>(receiver));
        return [receiverPointer, method](bool isEnteringState) {
            if (isEnteringState && !receiverPointer.isNull())
                (receiverPointer->*method)();
        };
    }

    template<typename PointerToMemberFunction>
    static std::function<void(bool)> onExit(
            const typename QtPrivate::FunctionPointer<PointerToMemberFunction>::Object *receiver,
            PointerToMemberFunction method)
    {
        typedef typename QtPrivate::FunctionPointer<PointerToMemberFunction>::Object Object;
        const QPointer<Object> receiverPointer(const_cast<Object *>(receiver));
        return [receiverPointer, method](bool isEnteringState) {
            if (!isEnteringState && !receiverPointer.isNull())
                (receiverPointer->*method)();
        };
    }

    QMetaObject::Connection connectToEvent(const QString &scxmlEventSpec,
                                           const QObject *receiver, const char *method,
                                           Qt::ConnectionType type = Qt::AutoConnection);

    // connect state to a QObject slot
    template <typename PointerToMemberFunction>
    inline QMetaObject::Connection connectToEvent(
            const QString &scxmlEventSpec,
            const typename QtPrivate::FunctionPointer<PointerToMemberFunction>::Object *receiver,
            PointerToMemberFunction method,
            Qt::ConnectionType type = Qt::AutoConnection)
    {
        typedef QtPrivate::FunctionPointer<PointerToMemberFunction> SlotType;
        return connectToEventImpl(
                    scxmlEventSpec, receiver, nullptr,
                    new QtPrivate::QSlotObject<PointerToMemberFunction,
                    typename SlotType::Arguments, void>(method),
                    type);
    }

    // connect state to a functor or function pointer (without context)
    template <typename Functor>
    inline typename QtPrivate::QEnableIf<
            !QtPrivate::FunctionPointer<Functor>::IsPointerToMemberFunction &&
            !std::is_same<const char*, Functor>::value, QMetaObject::Connection>::Type
    connectToEvent(const QString &scxmlEventSpec, Functor functor,
                   Qt::ConnectionType type = Qt::AutoConnection)
    {
        // Use this as context
        return connectToEvent(scxmlEventSpec, this, functor, type);
    }

    // connectToEvent to a functor or function pointer (with context)
    template <typename Functor>
    inline typename QtPrivate::QEnableIf<
            !QtPrivate::FunctionPointer<Functor>::IsPointerToMemberFunction &&
            !std::is_same<const char*, Functor>::value, QMetaObject::Connection>::Type
    connectToEvent(const QString &scxmlEventSpec, const QObject *context, Functor functor,
                   Qt::ConnectionType type = Qt::AutoConnection)
    {
        QtPrivate::QSlotObjectBase *slotObj = new QtPrivate::QFunctorSlotObject<Functor, 1,
                QtPrivate::List<QScxmlEvent>, void>(functor);
        return connectToEventImpl(scxmlEventSpec, context, reinterpret_cast<void **>(&functor),
                                  slotObj, type);
    }

    Q_INVOKABLE void submitEvent(QScxmlEvent *event);
    Q_INVOKABLE void submitEvent(const QString &eventName);
    Q_INVOKABLE void submitEvent(const QString &eventName, const QVariant &data);
    Q_INVOKABLE void cancelDelayedEvent(const QString &sendId);

    Q_INVOKABLE bool isDispatchableTarget(const QString &target) const;

    QVector<QScxmlInvokableService *> invokedServices() const;

    QScxmlTableData *tableData() const;
    void setTableData(QScxmlTableData *tableData);

Q_SIGNALS:
    void runningChanged(bool running);
    void invokedServicesChanged(const QVector<QScxmlInvokableService *> &invokedServices);
    void log(const QString &label, const QString &msg);
    void reachedStableState();
    void finished();
    void dataModelChanged(QScxmlDataModel *model);
    void initialValuesChanged(const QVariantMap &initialValues);
    void initializedChanged(bool initialized);
    void loaderChanged(QScxmlCompiler::Loader *loader);
    void tableDataChanged(QScxmlTableData *tableData);

public Q_SLOTS:
    void start();
    void stop();
    bool init();

protected: // methods for friends:
    friend class QScxmlDataModel;
    friend class QScxmlEventBuilder;
    friend class QScxmlInvokableServicePrivate;
    friend class QScxmlExecutionEngine;

    // The methods below are used by the compiled state machines.
    bool isActive(int stateIndex) const;

private:
    QMetaObject::Connection connectToStateImpl(const QString &scxmlStateName,
                                               const QObject *receiver, void **slot,
                                               QtPrivate::QSlotObjectBase *slotObj,
                                               Qt::ConnectionType type = Qt::AutoConnection);
    QMetaObject::Connection connectToEventImpl(const QString &scxmlEventSpec,
                                               const QObject *receiver, void **slot,
                                               QtPrivate::QSlotObjectBase *slotObj,
                                               Qt::ConnectionType type = Qt::AutoConnection);
};

QT_END_NAMESPACE

Q_DECLARE_METATYPE(QScxmlStateMachine *)
Q_DECLARE_METATYPE(QVector<QScxmlInvokableService *>)

#endif // QSCXMLSTATEMACHINE_H
