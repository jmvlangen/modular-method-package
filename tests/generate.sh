#!/bin/bash

if [ $# -ge 1 ]
then TEMPLATENAME=$1
else echo 'Please supply a template name:'
     read TEMPLATENAME
fi
if [[ "${TEMPLATENAME}" != *.sage ]]
then TEMPLATENAME="${TEMPLATENAME}.sage"
fi
TEMPLATEPATH="./templates/${TEMPLATENAME}"

if [ $# -ge 2 ]
then MIN=$2
else echo 'Please supply the minimal parameter:'
     read MIN
fi

if [ $# -ge 3 ]
then MAX=$3
else echo 'Please supply the maximal parameter:'
     read MAX
fi

VAR=$MIN
while [ $VAR -le $MAX ]
do  FILENAME="${TEMPLATENAME%.sage} ${VAR}.sage"
    FILEPATH="./code/${FILENAME}"
    cp "${TEMPLATEPATH}" "${FILEPATH}"
    replace '<n>' "${VAR}" -- "${FILEPATH}"
    ((VAR++))
done
