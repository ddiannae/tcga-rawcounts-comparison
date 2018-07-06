#!/bin/bash

LEGACYCANCER="legacy/cancer"
NEWCANCER="new/cancer"
LEGACYHEALTHY="legacy/healthy"
NEWHEALTHY="new/healthy"


find $NEWCANCER -name '*.gz' -exec mv '{}' $NEWCANCER \;
find $NEWHEALTHY -name '*.gz' -exec mv '{}' $NEWHEALTHY \;

find . -type d -empty -delete
find . -name '*.gz' -exec gunzip '{}' \;
find . -name '*.gz' -exec rm '{}' \;

find $LEGACYCANCER -name '*.gene.quantification.txt' -exec mv '{}' $LEGACYCANCER \;
find $LEGACYHEALTHY -name '*.gene.quantification.txt' -exec mv '{}' $LEGACYHEALTHY \;
