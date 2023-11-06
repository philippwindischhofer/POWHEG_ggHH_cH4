#!/bin/bash

(echo 1 ; ls -c1 pwg-*-NLO.top ; echo "") | $(dirname $0)/mergedata

mv fort.12 pwgNLO.top


