#!/bin/bash
set -euf -o pipefail

# Hardcoding this for the Scipy Dockerhub and scipy-dev image
REGISTRY=${REGISTRY:-"https://registry.hub.docker.com/v2/repositories/"}
IMAGE_NAME=${IMAGE_NAME:-"scipy/scipy-dev"}
TAGS_URL="${REGISTRY}${IMAGE_NAME}/tags"

# Some helpful vars
VERBOSE=1
ignore_404=0

echo "\xF0\x9F\x9A\xA9 Querying ${TAGS_URL}"

# Function to get the tags from the Docker registry
function do_curl_get() {
    local URL="$1"
    HTTP_RESPONSE="$(curl -sSL --write-out "HTTPSTATUS:%{http_code}" \
        -H "Content-Type: application/json;charset=UTF-8" \
        -X GET "$URL")"

    HTTP_BODY=$(echo "$HTTP_RESPONSE" | sed -E 's/HTTPSTATUS\:[0-9]{3}$//')
    HTTP_STATUS=$(echo "$HTTP_RESPONSE" | tr -d '\n' | sed -E 's/.*HTTPSTATUS:([0-9]{3})$/\1/')

    # Check that the http status is 200
    if [[ "$HTTP_STATUS" -ne 200 ]]; then
        if [[ "$ignore_404" -eq 0 ]]; then
            if [[ "$VERBOSE" -eq 0 ]]; then
                echo -e "\\nError $HTTP_STATUS from: $URL\\n"
            else
                echo -e "\\nError $HTTP_STATUS from: $URL\\nHTTP_BODY: $HTTP_BODY\\n"
            fi
        fi
    fi
}

# Get the tags from the public API
do_curl_get "$TAGS_URL"
TAGS_CURL=$(echo "$HTTP_BODY")

# Get the last updated tag - first remove latest from the tags and
# get the max (latest) push date -> extract the correstponding tag name
LAST_UPDATED=$(echo $TAGS_CURL | jq 'del(.results[] | select(.name == "latest")) |.results| max_by(.last_updated).name ' | sed 's/"//g')

echo $LAST_UPDATED

# Update the tag to the latest version
sed -E -i '' "s/[0-9]{8}\-master\-[a-z0-9]{8}/${LAST_UPDATED}/" tools/docker_dev/.gitpod.Dockerfile
