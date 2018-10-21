#!/bin/bash

#for i in $(seq -f%04g 6 76)
for i in $(seq -f%04g 77 93)
do
#  echo "./r2root /run/media/jing.li/May312018/user/1769_data1x/Run00${i}HFC/HFC.dat ../1769_root/skim00${i}.root CC.cal ../1769_skim/skim00${i}.dat"
  #./r2root /run/media/jing.li/June222018/user/1563_track/TK${i}.dat ../1563_root/run0${i}.root CC.cal
  ./r2root ../1563_track/TK${i}.dat ../1563_root/run${i}.root CC.cal
done
