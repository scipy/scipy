/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtBluetooth module of the Qt Toolkit.
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

#ifndef QPRIVATELINEARBUFFER_P_H
#define QPRIVATELINEARBUFFER_P_H

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

// This is QIODevice's read buffer, optimised for read(), isEmpty() and getChar()
class QPrivateLinearBuffer
{
public:
    QPrivateLinearBuffer() : len(0), first(0), buf(0), capacity(0) {
    }
    ~QPrivateLinearBuffer() {
        delete [] buf;
    }
    void clear() {
        first = buf;
        len = 0;
    }
    int size() const {
        return len;
    }
    bool isEmpty() const {
        return len == 0;
    }
    void skip(int n) {
        if (n >= len) {
            clear();
        } else {
            len -= n;
            first += n;
        }
    }
    int getChar() {
        if (len == 0)
            return -1;
        int ch = uchar(*first);
        len--;
        first++;
        return ch;
    }
    int read(char* target, int size) {
        int r = qMin(size, len);
        memcpy(target, first, r);
        len -= r;
        first += r;
        return r;
    }
    char* reserve(int size) {
        makeSpace(size + len, freeSpaceAtEnd);
        char* writePtr = first + len;
        len += size;
        return writePtr;
    }
    void chop(int size) {
        if (size >= len) {
            clear();
        } else {
            len -= size;
        }
    }
    QByteArray readAll() {
        char* f = first;
        int l = len;
        clear();
        return QByteArray(f, l);
    }
    int readLine(char* target, int size) {
        int r = qMin(size, len);
        char* eol = static_cast<char*>(memchr(first, '\n', r));
        if (eol)
            r = 1+(eol-first);
        memcpy(target, first, r);
        len -= r;
        first += r;
        return int(r);
        }
    bool canReadLine() const {
        return memchr(first, '\n', len);
    }
    void ungetChar(char c) {
        if (first == buf) {
            // underflow, the existing valid data needs to move to the end of the (potentially bigger) buffer
            makeSpace(len+1, freeSpaceAtStart);
        }
        first--;
        len++;
        *first = c;
    }
    void ungetBlock(const char* block, int size) {
        if ((first - buf) < size) {
            // underflow, the existing valid data needs to move to the end of the (potentially bigger) buffer
            makeSpace(len + size, freeSpaceAtStart);
        }
        first -= size;
        len += size;
        memcpy(first, block, size);
    }

private:
    enum FreeSpacePos {freeSpaceAtStart, freeSpaceAtEnd};
    void makeSpace(size_t required, FreeSpacePos where) {
        size_t newCapacity = qMax(capacity, size_t(QPRIVATELINEARBUFFER_BUFFERSIZE));
        while (newCapacity < required)
            newCapacity *= 2;
        const int moveOffset = (where == freeSpaceAtEnd) ? 0 : int(newCapacity) - len;
        if (newCapacity > capacity) {
            // allocate more space
            char* newBuf = new char[newCapacity];
            memmove(newBuf + moveOffset, first, len);
            delete [] buf;
            buf = newBuf;
            capacity = newCapacity;
        } else {
            // shift any existing data to make space
            memmove(buf + moveOffset, first, len);
        }
        first = buf + moveOffset;
    }

private:
    // length of the unread data
    int len;
    // start of the unread data
    char* first;
    // the allocated buffer
    char* buf;
    // allocated buffer size
    size_t capacity;
};

#endif // QPRIVATELINEARBUFFER_P_H
