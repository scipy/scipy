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

#ifndef QSCXMLSTATEMACHINE_P_H
#define QSCXMLSTATEMACHINE_P_H

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

#include <QtScxml/private/qscxmlexecutablecontent_p.h>
#include <QtScxml/qscxmlstatemachine.h>
#include <QtScxml/private/qscxmlstatemachineinfo_p.h>
#include <QtCore/private/qobject_p.h>
#include <QtCore/private/qmetaobject_p.h>
#include <QtCore/qmetaobject.h>

QT_BEGIN_NAMESPACE

namespace QScxmlInternal {
class EventLoopHook: public QObject
{
    Q_OBJECT

    QScxmlStateMachinePrivate *smp;

public:
    EventLoopHook(QScxmlStateMachinePrivate *smp)
        : smp(smp)
    {}

    void queueProcessEvents();

    Q_INVOKABLE void doProcessEvents();

protected:
    void timerEvent(QTimerEvent *timerEvent) override;
};

class ScxmlEventRouter : public QObject
{
    Q_OBJECT
public:
    ScxmlEventRouter(QObject *parent = nullptr) : QObject(parent) {}
    QMetaObject::Connection connectToEvent(const QStringList &segments, const QObject *receiver,
                                           const char *method, Qt::ConnectionType type);
    QMetaObject::Connection connectToEvent(const QStringList &segments, const QObject *receiver,
                                           void **slot, QtPrivate::QSlotObjectBase *method,
                                           Qt::ConnectionType type);

    void route(const QStringList &segments, QScxmlEvent *event);

signals:
    void eventOccurred(const QScxmlEvent &event);

private:
    QHash<QString, ScxmlEventRouter *> children;
    ScxmlEventRouter *child(const QString &segment);

    void disconnectNotify(const QMetaMethod &signal) override;
};

class StateMachineInfoProxy: public QObject
{
    Q_OBJECT

public:
    StateMachineInfoProxy(QObject *parent)
        : QObject(parent)
    {}

Q_SIGNALS:
    void statesEntered(const QVector<QScxmlStateMachineInfo::StateId> &states);
    void statesExited(const QVector<QScxmlStateMachineInfo::StateId> &states);
    void transitionsTriggered(const QVector<QScxmlStateMachineInfo::TransitionId> &transitions);
};
} // QScxmlInternal namespace

class QScxmlInvokableService;
class QScxmlStateMachinePrivate: public QObjectPrivate
{
    Q_DECLARE_PUBLIC(QScxmlStateMachine)

    static QAtomicInt m_sessionIdCounter;

public: // types
    typedef QScxmlExecutableContent::StateTable StateTable;

    class HistoryContent
    {
        QHash<int, int> storage;

    public:
        HistoryContent() { storage.reserve(4); }

        int &operator[](int idx) {
            QHash<int, int>::Iterator i = storage.find(idx);
            return (i == storage.end()) ? storage.insert(idx, StateTable::InvalidIndex).value() :
                                          i.value();
        }

        int value(int idx) const {
            QHash<int, int>::ConstIterator i = storage.constFind(idx);
            return (i == storage.constEnd()) ? StateTable::InvalidIndex : i.value();
        }
    };

    class ParserData
    {
    public:
        QScopedPointer<QScxmlDataModel> m_ownedDataModel;
        QVector<QScxmlError> m_errors;
    };

    // The OrderedSet is a set where it elements are in insertion order. See
    // http://www.w3.org/TR/scxml/#AlgorithmforSCXMLInterpretation under Algorithm, Datatypes. It
    // is used to keep lists of states and transitions in document order.
    class OrderedSet
    {
        std::vector<int> storage;

    public:
        OrderedSet(){}
        OrderedSet(std::initializer_list<int> l): storage(l) {}

        std::vector<int> takeList() const
        { return std::move(storage); }

        const std::vector<int> &list() const
        { return storage; }

        bool contains(int i) const
        {
            return std::find(storage.cbegin(), storage.cend(), i) != storage.cend();
        }

        bool remove(int i)
        {
            std::vector<int>::iterator it = std::find(storage.begin(), storage.end(), i);
            if (it == storage.end()) {
                return false;
            }
            storage.erase(it);
            return true;
        }

        void removeHead()
        { if (!isEmpty()) storage.erase(storage.begin()); }

        bool isEmpty() const
        { return storage.empty(); }

        void add(int i)
        { if (!contains(i)) storage.push_back(i); }

        bool intersectsWith(const OrderedSet &other) const
        {
            for (auto i : storage) {
                if (other.contains(i)) {
                    return true;
                }
            }
            return false;
        }

        void clear()
        { storage.clear(); }

        typedef std::vector<int>::const_iterator const_iterator;
        const_iterator begin() const { return storage.cbegin(); }
        const_iterator end() const { return storage.cend(); }
    };

    class Queue
    {
        QVector<QScxmlEvent *> storage;

    public:
        Queue()
        { storage.reserve(4); }

        ~Queue()
        {  qDeleteAll(storage); }

        void enqueue(QScxmlEvent *e)
        { storage.append(e); }

        bool isEmpty() const
        { return storage.empty(); }

        QScxmlEvent *dequeue()
        {
            Q_ASSERT(!isEmpty());
            QScxmlEvent *e = storage.first();
            storage.pop_front();
            int sz = storage.size();
            if (Q_UNLIKELY(sz > 4 && sz * 8 < storage.capacity())) {
                storage.squeeze();
            }
            return e;
        }
    };

public:
    QScxmlStateMachinePrivate(const QMetaObject *qMetaObject);
    ~QScxmlStateMachinePrivate();

    static QScxmlStateMachinePrivate *get(QScxmlStateMachine *t)
    { return t->d_func(); }

    static QString generateSessionId(const QString &prefix);

    ParserData *parserData();

    void setIsInvoked(bool invoked)
    { m_isInvoked = invoked; }

    void addService(int invokingState);
    void removeService(int invokingState);
    QScxmlInvokableServiceFactory *serviceFactory(int id);

    bool executeInitialSetup();

    void routeEvent(QScxmlEvent *event);
    void postEvent(QScxmlEvent *event);
    void submitDelayedEvent(QScxmlEvent *event);
    void submitError(const QString &type, const QString &msg, const QString &sendid = QString());

    void start();
    void pause();
    void processEvents();

    void setEvent(QScxmlEvent *event);
    void resetEvent();

    void emitStateActive(int stateIndex, bool active);
    void emitInvokedServicesChanged();

    void attach(QScxmlStateMachineInfo *info);
    const OrderedSet &configuration() const { return m_configuration; }

    void updateMetaCache();

private:
    QStringList stateNames(const std::vector<int> &stateIndexes) const;
    std::vector<int> historyStates(int stateIdx) const;

    void exitInterpreter();
    void returnDoneEvent(QScxmlExecutableContent::ContainerId doneData);
    bool nameMatch(const StateTable::Array &patterns, QScxmlEvent *event) const;
    void selectTransitions(OrderedSet &enabledTransitions,
                           const std::vector<int> &configInDocumentOrder,
                           QScxmlEvent *event) const;
    void removeConflictingTransitions(OrderedSet *enabledTransitions) const;
    void getProperAncestors(std::vector<int> *ancestors, int state1, int state2) const;
    void microstep(const OrderedSet &enabledTransitions);
    void exitStates(const OrderedSet &enabledTransitions);
    void computeExitSet(const OrderedSet &enabledTransitions, OrderedSet &statesToExit) const;
    void executeTransitionContent(const OrderedSet &enabledTransitions);
    void enterStates(const OrderedSet &enabledTransitions);
    void computeEntrySet(const OrderedSet &enabledTransitions,
                         OrderedSet *statesToEnter,
                         OrderedSet *statesForDefaultEntry,
                         HistoryContent *defaultHistoryContent) const;
    void addDescendantStatesToEnter(int stateIndex,
                                    OrderedSet *statesToEnter,
                                    OrderedSet *statesForDefaultEntry,
                                    HistoryContent *defaultHistoryContent) const;
    void addAncestorStatesToEnter(int stateIndex,
                                  int ancestorIndex,
                                  OrderedSet *statesToEnter,
                                  OrderedSet *statesForDefaultEntry,
                                  HistoryContent *defaultHistoryContent) const;
    std::vector<int> getChildStates(const StateTable::State &state) const;
    bool hasDescendant(const OrderedSet &statesToEnter, int childIdx) const;
    bool allDescendants(const OrderedSet &statesToEnter, int childdx) const;
    bool isDescendant(int state1, int state2) const;
    bool allInFinalStates(const std::vector<int> &states) const;
    bool someInFinalStates(const std::vector<int> &states) const;
    bool isInFinalState(int stateIndex) const;
    int getTransitionDomain(int transitionIndex) const;
    int findLCCA(OrderedSet &&states) const;
    void getEffectiveTargetStates(OrderedSet *targets, int transitionIndex) const;

public: // types & data fields:
    QString m_sessionId;
    bool m_isInvoked;
    bool m_isInitialized;
    bool m_isProcessingEvents;
    QVariantMap m_initialValues;
    QScxmlDataModel *m_dataModel;
    QScxmlCompilerPrivate::DefaultLoader m_defaultLoader;
    QScxmlCompiler::Loader *m_loader;
    QScxmlExecutionEngine *m_executionEngine;
    QScxmlTableData *m_tableData;
    const StateTable *m_stateTable;
    QScxmlStateMachine *m_parentStateMachine;
    QScxmlInternal::EventLoopHook m_eventLoopHook;
    typedef std::vector<std::pair<int, QScxmlEvent *>> DelayedQueue;
    DelayedQueue m_delayedEvents;
    const QMetaObject *m_metaObject;
    QScxmlInternal::ScxmlEventRouter m_router;

private:
    QScopedPointer<ParserData> m_parserData; // used when created by StateMachine::fromFile.
    typedef QHash<int, QVector<int>> HistoryValues;
    struct InvokedService {
        int invokingState;
        QScxmlInvokableService *service;
        QString serviceName;
    };

    // TODO: move the stuff below to a struct that can be reset
    HistoryValues m_historyValue;
    OrderedSet m_configuration;
    Queue m_internalQueue;
    Queue m_externalQueue;
    QSet<int> m_statesToInvoke;
    std::vector<InvokedService> m_invokedServices;
    std::vector<bool> m_isFirstStateEntry;
    std::vector<QScxmlInvokableServiceFactory *> m_cachedFactories;
    enum { Invalid = 0, Starting, Running, Paused, Finished } m_runningState = Invalid;
    bool isRunnable() const {
        switch (m_runningState) {
        case Starting:
        case Running:
        case Paused:
            return true;
        case Invalid:
        case Finished:
            return false;
        }

        return false; // Dead code, but many dumb compilers cannot (or are unwilling to) detect that.
    }

    bool isPaused() const { return m_runningState == Paused; }

    QScxmlInternal::StateMachineInfoProxy *m_infoSignalProxy;

    QHash<int, int> m_stateIndexToSignalIndex;
    QHash<QString, int> m_stateNameToSignalIndex;
};

QT_END_NAMESPACE

#endif // QSCXMLSTATEMACHINE_P_H

