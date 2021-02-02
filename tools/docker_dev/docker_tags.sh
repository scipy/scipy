#!/bin/bash

#!/bin/bash

if [ $# -lt 1 ]
then
cat << HELP

dockertags  --  list all tags for a Docker image on a remote registry.

EXAMPLE: 
    - list all tags for scipy:
       dockertags scipy-dev

    - list all scipy tags containing dev:
       dockertags scipy dev

HELP
fi

image="$1"
tags=`http https://registry.hub.docker.com/v1/repositories/scipy/${image}/tags   | sed -e 's/[][]//g' -e 's/"//g' -e 's/ //g' | tr '}' '\n'  | awk -F: '{print $3}'`


curl -L -s 'https://registry.hub.docker.com/v2/repositories/library/centos/tags?page_size=1024'|jq '."results"[]["name"]'

http https://registry.hub.docker.com/v2/repositories/scipy/scipy-dev/tags | jq '.results[] | [.name, .tag_last_pushed]'
http https://registry.hub.docker.com/v2/repositories/scipy/scipy-dev/tags | jq '.results[] | {tag: .name, date: .tag_last_pushed }'

http https://registry.hub.docker.com/v2/repositories/scipy/scipy-dev/tags | jq '.results |= sort_by(.last_updated)'
http https://registry.hub.docker.com/v2/repositories/scipy/scipy-dev/tags | jq '.results |= max_by(.last_updated)'

http https://registry.hub.docker.com/v2/repositories/scipy/scipy-dev/tags | jq '.results |= max_by(.last_updated) | .results | .name'

if [ -n "$2" ]
then
    tags=` echo "${tags}" | grep "$2" `
fi

echo "${tags}"