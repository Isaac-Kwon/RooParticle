#!/bin/bash

curl --request POST --header "PRIVATE-TOKEN:$GITLAB_API_KEY" --form "note=Nice picture man!" https://gitlab.cern.ch/api/v4/projects/${CI_PROJECT_PATH//\//%2F}/repository/commits/$CI_BUILD_REF/comments