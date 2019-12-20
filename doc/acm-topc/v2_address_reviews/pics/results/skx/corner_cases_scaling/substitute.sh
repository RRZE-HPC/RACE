#!/bin/bash

#USAGE: ./substitute.sh old_string new_string file

old_str="%$1%"
new_str="$2"
file=$3

sed -i "s#$old_str#$new_str#g" $file

