#!/bin/bash

#for i in $(seq -f%04g 1 76)
for i in $(seq -f%04g 20 57)
do
#  echo "./r2root /run/media/jing.li/May312018/user/1769_data1x/Run00${i}HFC/HFC.dat ../1769_root/skim00${i}.root CC.cal ../1769_skim/skim00${i}.dat"
  #./r2root /run/media/jing.li/June222018/user/1563_track/TK${i}.dat ../1563_root/run${i}.root CC.cal
  ./r2root /run/media/jing.li/June142018/user/1769_track/TK${i}.dat ../1769_root/tk${i}.root CC2.cal
done
