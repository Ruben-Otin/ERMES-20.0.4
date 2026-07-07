#!/bin/bash

# ProblemName     = $1
# ProblemPath     = $2
# ProblemTypePath = $3

# OutputFile: $2/$1.info

rm -f "$2/$1.info" 
rm -f "$2/$1.flavia.res" 
rm -f "$2/$1.post.res"
rm -f "$2/$1.post.msh"

"$3/execs/unix/ERMES_20.0.4_c7.5" "$1" > "$2/$1.info"