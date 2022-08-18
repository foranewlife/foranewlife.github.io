#!/bin/zsh
for i in ${1}/*
do
  python encode.py $i
done