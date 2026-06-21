#!/bin/bash

# SPDX-FileCopyrightText: 2020 Intel Corporation
#
# SPDX-License-Identifier: MIT

URL=$1
COMPONENTS=$2
EXPECTED_SHA256=$3

curl --output webimage.sh --url "$URL" --retry 5 --retry-delay 5

if [ -n "$EXPECTED_SHA256" ]; then
  file_hash=$(sha256sum webimage.sh | cut -d ' ' -f 1)
  if [ "$file_hash" != "$EXPECTED_SHA256" ]; then
    echo "Checksum verification failed for $URL"
    echo "Expected: $EXPECTED_SHA256"
    echo "Got:      $file_hash"
    exit 1
  fi
fi

chmod +x webimage.sh
./webimage.sh -x -f webimage_extracted --log extract.log
rm -rf webimage.sh
WEBIMAGE_NAME=$(ls -1 webimage_extracted/)
if [ -z "$COMPONENTS" ]; then
  sudo webimage_extracted/"$WEBIMAGE_NAME"/bootstrapper -s --action install --eula=accept --log-dir=.
  installer_exit_code=$?
else
  sudo webimage_extracted/"$WEBIMAGE_NAME"/bootstrapper -s --action install --components="$COMPONENTS" --eula=accept --log-dir=.
  installer_exit_code=$?
fi
rm -rf webimage_extracted
exit $installer_exit_code
