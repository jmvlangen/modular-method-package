#!/bin/bash
# Cleans up the results directory and the code directory if -c is set.

if [ "$1" == '-c' ]
then rm ./code/*
fi
rm ./results/*
