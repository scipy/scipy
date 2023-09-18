# Copyright (c) 2022, Riverbank Computing Limited
# All rights reserved.
#
# This copy of SIP is licensed for use under the terms of the SIP License
# Agreement.  See the file LICENSE for more details.
#
# This copy of SIP may also used under the terms of the GNU General Public
# License v2 or v3 as published by the Free Software Foundation which can be
# found in the files LICENSE-GPL2 and LICENSE-GPL3 included in this package.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ('AS IS'
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.


from ...exceptions import UserException


def output_extract(spec, extract_ref):
    """ Output an extract. """

    # Get the id and file name from the reference.
    extract_id, _, extract_file = extract_ref.partition(':')
    if not extract_file:
        raise UserException(
                "an extract must be in the form 'id:file', not '{0}'".format(
                        extract_ref))

    # Get all the parts of this extract.
    ordered_parts = []
    appended_parts = []

    for extract in spec.extracts:
        if extract.id == extract_id:
            if extract.order < 0:
                appended_parts.append(extract)
            else:
                ordered_parts.append(extract)

    # Put the parts in the correct order.
    ordered_parts.sort(key=lambda e: e.order)
    ordered_parts.extend(appended_parts)

    if len(ordered_parts) == 0:
        raise UserException(
                "there is no extract defined with the identifier '{0}'".format(
                        extract_id))

    # Write the parts to the file.
    with open(extract_file, 'w') as ef:
        for extract_part in ordered_parts:
            ef.write(extract_part.text)
