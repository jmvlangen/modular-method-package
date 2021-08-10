#!/usr/bin/env bash

# The help function
show_help() {
    echo "Compiles all sage source files into a python module"
    echo ""
    echo "Usage:"
    echo " >> ./compile.sh"
    echo ""
    echo "Options:"
    echo " -f force replace all files"
    echo " -i interactively replace all files"
    echo " -n do not replace any files"
    echo " --help show this text and exit"
    exit 0
}

# Parsing options
replace=1 # Replace older files only
while [[ -n $1 ]]; do
    case $1 in
	-f) replace=3;; # Force replace
	-i) replace=2;; # Interactive replace
	-n) replace=0;; # No replacing
	--help) show_help;;
	*) echo "Option $1 is not recognized";;
    esac
    shift
done

ask_replace() {
    echo "Do you wish to replace $1? Y(es) or N(o): "
    while read answer; do
	case $answer in
	    Y|y|Yes|yes) return 0;;
	    N|n|No|no) return 1;;
	    *) echo "Replace $1? Please answer Y(es) or N(o): ";;
	esac
    done
}

compile_file() {
    name=${1%.sage}
    if [[ ! -e "${name}.py" ]] ||
	   [[ $replace == 3 ]] ||
	   { [[ $replace == 2 ]] && ask_replace "${name}.py"; } ||
	   { [[ $replace == 1 ]] && [[ "${name}.sage" -nt "${name}.py" ]]; }; then
	if [[ -e "${name}.py" ]]; then
	    rm "${name}.py"
	fi
	echo "Compiling ${name}.sage"
	sage -preparse "${name}.sage"
	mv "${name}.sage.py" "${name}.py"
	sed -i -f commands.sed "${name}.py"
    fi
    rm "${name}.sage"
}

compile_directory() {
    if [[ "$1" =~ .*/__pycache__ ]]; then
	return
    fi
    echo "Compiling ${1}"
    for file in ${1}/*; do
	if [[ -d $file ]]; then
	    compile_directory "$file"
	elif [[ "$file" =~ .+\.sage ]]; then
	    compile_file "$file"
	fi
    done
    if ! [[ -e "${1}/__init__.py" ]]; then
	touch "${1}/__init__.py"
    fi
}

target="modular_method"
if ! [[ -e "$target" ]] || ! [[ -d "$target" ]]; then
    mkdir "$target"
fi
cp -r --preserve=timestamps src/* modular_method/
compile_directory "$target"
