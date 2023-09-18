// Copyright 2020 The Chromium Authors. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file.

/**
 * @fileoverview Rule to check license headers
 * @author Tim van der Lippe
 */
'use strict';

const path = require('path');

const FRONT_END_FOLDER = path.join(__filename, '..', '..', '..', '..', 'front_end');

const LINE_LICENSE_HEADER = [
  'Copyright 2020 The Chromium Authors. All rights reserved.',
  'Use of this source code is governed by a BSD-style license that can be',
  'found in the LICENSE file.',
];

const BLOCK_LICENSE_HEADER = [
  'Copyright \\(C\\) \\d{4} Google Inc. All rights reserved.',
  '',
  'Redistribution and use in source and binary forms, with or without',
  'modification, are permitted provided that the following conditions are',
  'met:',
  '',
  '    \\* Redistributions of source code must retain the above copyright',
  'notice, this list of conditions and the following disclaimer.',
  '    \\* Redistributions in binary form must reproduce the above',
  'copyright notice, this list of conditions and the following disclaimer',
  'in the documentation and/or other materials provided with the',
  'distribution.',
  '    \\* Neither the name of Google Inc. nor the names of its',
  'contributors may be used to endorse or promote products derived from',
  'this software without specific prior written permission.',
  '',
  'THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS',
  '"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT',
  'LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR',
  'A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT',
  'OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,',
  'SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES \\(INCLUDING, BUT NOT',
  'LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,',
  'DATA, OR PROFITS; OR BUSINESS INTERRUPTION\\) HOWEVER CAUSED AND ON ANY',
  'THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT',
  '\\(INCLUDING NEGLIGENCE OR OTHERWISE\\) ARISING IN ANY WAY OUT OF THE USE',
  'OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.',
];

const LINE_REGEXES = LINE_LICENSE_HEADER.map(line => new RegExp('[ ]?' + line.replace('2020', '(\\(c\\) )?\\d{4}')));
const BLOCK_REGEX = new RegExp('[\\s\\\\n\\*]*' + BLOCK_LICENSE_HEADER.join('[\\s\\\\n\\*]*'), 'm');

const LICENSE_HEADER_ADDITION = LINE_LICENSE_HEADER.map(line => `// ${line}`).join('\n') + '\n\n';

const EXCLUDED_FILES = [
  // FIXME: Dagre bundles must be moved to third_party
  'dagre_layout/dagre.js',
  // FIXME: Diff bundles must be moved to third_party
  'diff/diff_match_patch.js',
];

const OTHER_LICENSE_HEADERS = [
  // Apple
  'bindings/ResourceUtils.js',
  'common/Color.js',
  'common/Object.js',
  'common/ResourceType.js',
  'data_grid/DataGrid.js',
  'dom_extension/DOMExtension.js',
  'elements/MetricsSidebarPane.js',
  'profiler/CPUProfileView.js',
  'profiler/ProfilesPanel.js',
  'resources/ApplicationCacheItemsView.js',
  'resources/ApplicationCacheModel.js',
  'resources/DatabaseModel.js',
  'resources/DatabaseQueryView.js',
  'resources/DatabaseTableView.js',
  'sdk/Resource.js',
  'sdk/Script.js',
  'source_frame/FontView.js',
  'source_frame/ImageView.js',
  'sources/CallStackSidebarPane.js',
  'ui/Panel.js',
  'ui/Treeoutline.js',
  // Brian Grinstead
  'color_picker/Spectrum.js',
  // Joseph Pecoraro
  'console/ConsolePanel.js',
  // Research In Motion Limited
  'network/ResourceWebSocketFrameView.js',
  // 280 North Inc.
  'profiler/BottomUpProfileDataGrid.js',
  'profiler/ProfileDataGrid.js',
  'profiler/TopDownProfileDataGrid.js',
  // IBM Corp
  'sources/WatchExpressionsSidebarPane.js',
  // Multiple authors
  'components/JSPresentationUtils.js',
  'console/ConsoleView.js',
  'console/ConsoleViewMessage.js',
  'cookie_table/CookiesTable.js',
  'elements/ComputedStyleWidget.js',
  'elements/ElementsPanel.js',
  'elements/ElementsTreeElement.js',
  'elements/ElementsTreeOutline.js',
  'elements/EventListenersWidget.js',
  'elements/PropertiesWidget.js',
  'elements/StylesSidebarPane.js',
  'main/MainImpl.js',
  'network/HARWriter.js',
  'network/NetworkDataGridNode.js',
  'network/NetworkLogView.js',
  'network/NetworkPanel.js',
  'network/NetworkTimeCalculator.js',
  'network/RequestHeadersView.js',
  'object_ui/ObjectPropertiesSection.js',
  'perf_ui/TimelineGrid.js',
  'platform/utilities.js',
  'platform/UIString.js',
  'resources/ApplicationPanelSidebar.js',
  'resources/CookieItemsView.js',
  'resources/DOMStorageItemsView.js',
  'resources/DOMStorageModel.js',
  'sdk/DOMModel.js',
  'source_frame/ResourceSourceFrame.js',
  'sources/ScopeChainSidebarPane.js',
  'sources/SourcesPanel.js',
  'theme_support/theme_support_impl.js',
  'timeline/TimelinePanel.js',
  'timeline/TimelineUIUtils.js',
  'ui/KeyboardShortcut.js',
  'ui/SearchableView.js',
  'ui/TextPrompt.js',
  'ui/UIUtils.js',
  'ui/Widget.js',
];

// ------------------------------------------------------------------------------
// Rule Definition
// ------------------------------------------------------------------------------

/**
 * Check each linecomment that should (combined) result in the LINE_LICENSE_HEADER.
 */
function isMissingLineCommentLicense(comments) {
  for (let i = 0; i < LINE_REGEXES.length; i++) {
    if (!comments[i] || !LINE_REGEXES[i].test(comments[i].value)) {
      return true;
    }
  }

  return false;
}

/**
 * We match the whole block comment, including potential leading asterisks of the jsdoc.
 */
function isMissingBlockLineCommentLicense(licenseText) {
  return !BLOCK_REGEX.test(licenseText);
}

module.exports = {
  meta: {
    type: 'problem',

    docs: {
      description: 'check license headers',
      category: 'Possible Errors',
    },
    fixable: 'code',
    schema: []  // no options
  },
  create: function(context) {
    const fileName = context.getFilename();
    // Fix windows paths for exemptions
    const relativePath = path.relative(FRONT_END_FOLDER, fileName).replace(/\\/g, '/');

    if (relativePath.startsWith('third_party') || fileName.endsWith('TestRunner.js') ||
        EXCLUDED_FILES.includes(relativePath) || OTHER_LICENSE_HEADERS.includes(relativePath)) {
      return {};
    }

    return {
      Program(node) {
        if (node.body.length === 0) {
          return;
        }

        const {leading: comments} = context.getComments(node.body[0]);

        if (!comments || comments.length === 0) {
          context.report({
            node,
            message: 'Missing license header',
            fix(fixer) {
              return fixer.insertTextBefore(node, LICENSE_HEADER_ADDITION);
            },
          });
        } else if (comments[0].type === 'Line') {
          if (isMissingLineCommentLicense(comments)) {
            context.report({
              node,
              message: 'Incorrect line license header',
              fix(fixer) {
                return fixer.insertTextBefore(comments[0], LICENSE_HEADER_ADDITION);
              }
            });
          }
        } else {
          if (isMissingBlockLineCommentLicense(comments[0].value)) {
            context.report({
              node,
              message: 'Incorrect block license header',
              fix(fixer) {
                return fixer.insertTextBefore(comments[0], LICENSE_HEADER_ADDITION);
              }
            });
          }
        }
      }
    };
  }
};
