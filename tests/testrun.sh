#!/bin/bash
timetmp="./results/time.tmp"
touch "${timetmp}"

for path in ./code/*.sage
do
    file=${path##*/}
    name=${file%.sage}
    output="./results/${name}.log"
    touch "${output}"
    /usr/bin/time -v -o "${timetmp}" sage "${path}" &> "${output}"
    echo '' >> "${output}"
    echo '--------------------------------------------------------------' >> "${output}"
    echo '' >> "${output}"
    cat "${timetmp}" >> "${output}"
done

rm "${timetmp}"
