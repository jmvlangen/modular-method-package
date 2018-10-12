#!/bin/bash
timetmp="./results/time.tmp"
runsource="./templates/run.sage"
runtmp="./results/run.sage"
logtmp="./results/run.log"
touch "${runtmp}"
touch "${timetmp}"
touch "${logtmp}"

IFS=''
find ~/Documents/SageFiles/database -name '*.sage' -printf "%p\n" | while read path;
do
    cp "${runsource}" "${runtmp}"
    replace '<path>' "${path}" -- "${runtmp}"
    echo "Running ${path}"
    touch "${logtmp}"
    /usr/bin/time -v -o "${timetmp}" sage "${runtmp}" &> "${logtmp}"
    echo '' >> "${logtmp}"
    echo '--------------------------------------------------------------' >> "${logtmp}"
    echo '' >> "${logtmp}"
    cat "${timetmp}" >> "${logtmp}"
    mv -f "${logtmp}" "${path%.sage}.log"
done

rm "${runtmp}"
rm "${runtmp}.py"
rm "${logtmp}"
rm "${timetmp}"
