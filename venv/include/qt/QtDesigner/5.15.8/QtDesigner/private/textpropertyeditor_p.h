/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt Designer of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:GPL-EXCEPT$
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
** General Public License version 3 as published by the Free Software
** Foundation with exceptions as appearing in the file LICENSE.GPL3-EXCEPT
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
// This file is not part of the Qt API.  It exists for the convenience
// of Qt Designer.  This header
// file may change from version to version without notice, or even be removed.
//
// We mean it.
//

#ifndef TEXTPROPERTYEDITOR_H
#define TEXTPROPERTYEDITOR_H

#include "shared_global_p.h"
#include "shared_enums_p.h"

#include <QtWidgets/qwidget.h>

QT_BEGIN_NAMESPACE

namespace qdesigner_internal {

    class PropertyLineEdit;

    // Inline-Editor for text properties. Does escaping of newline characters
    // to '\n' and back and provides validation modes. The interface
    // corresponds to that of QLineEdit.
    class QDESIGNER_SHARED_EXPORT TextPropertyEditor : public QWidget
    {
        TextPropertyEditor(const TextPropertyEditor &);
        TextPropertyEditor& operator=(const TextPropertyEditor &);
        Q_OBJECT
        Q_PROPERTY(QString text READ text WRITE setText USER true)
    public:
        enum EmbeddingMode {
            // Stand-alone widget
            EmbeddingNone,
                // Disable frame
                EmbeddingTreeView,
                // For editing in forms
                EmbeddingInPlace
        };

        enum UpdateMode {
            // Emit textChanged() as the user types
            UpdateAsYouType,
            // Emit textChanged() only when the user finishes (for QUrl, etc.)
            UpdateOnFinished
        };

        explicit TextPropertyEditor(QWidget *parent = nullptr, EmbeddingMode embeddingMode = EmbeddingNone, TextPropertyValidationMode validationMode = ValidationMultiLine);

        TextPropertyValidationMode textPropertyValidationMode() const { return m_validationMode; }
        void setTextPropertyValidationMode(TextPropertyValidationMode vm);

        UpdateMode updateMode() const                { return m_updateMode; }
        void setUpdateMode(UpdateMode um) { m_updateMode = um; }

        QString text() const;

        QSize sizeHint () const override;
        QSize minimumSizeHint () const override;

        void setAlignment(Qt::Alignment alignment);

        bool hasAcceptableInput() const;

        // installs an event filter object on the private QLineEdit
        void installEventFilter(QObject *filterObject);

        // Replace newline characters by literal "\n" for inline editing
        // in mode ValidationMultiLine
        static QString stringToEditorString(const QString &s, TextPropertyValidationMode validationMode = ValidationMultiLine);

        // Replace literal "\n"  by actual new lines in mode ValidationMultiLine
        static QString editorStringToString(const QString &s, TextPropertyValidationMode validationMode = ValidationMultiLine);

        // Returns whether newline characters are valid in validationMode.
        static bool multiLine(TextPropertyValidationMode validationMode);

    signals:
        void textChanged(const QString &text);
        void editingFinished();

    public slots:
        void setText(const QString &text);
        void selectAll();
        void clear();

    protected:
        void resizeEvent(QResizeEvent * event ) override;

    private slots:
        void slotTextChanged(const QString &text);
        void slotTextEdited();
        void slotEditingFinished();

    private:
        void setRegularExpressionValidator(const QString &pattern);
        void markIntermediateState();

        TextPropertyValidationMode m_validationMode = ValidationSingleLine;
        UpdateMode m_updateMode = UpdateAsYouType;
        PropertyLineEdit* m_lineEdit;

        // Cached text containing real newline characters.
        QString m_cachedText;
        bool m_textEdited = false;
    };
}

QT_END_NAMESPACE

#endif // TEXTPROPERTYEDITOR_H
