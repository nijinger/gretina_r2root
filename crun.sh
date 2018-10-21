#!/bin/bash

for i in $(seq 177 200)
do
#  echo "./r2root /run/media/jing.li/May312018/user/1769_data1x/Run00${i}HFC/HFC.dat ../1769_root/skim00${i}.root CC.cal ../1769_skim/skim00${i}.dat"
#  ./r2root /run/media/jing.li/May312018/user/1769_data1x/Run00${i}HFC/HFC.dat ../1769_root/skim00${i}.root CC.cal ../1769_skim/skim00${i}.dat
#  ./r2root /run/media/jing.li/June162018/user/1769_data1x/Run0${i}HFC/HFC.dat ../1769_root/skim0${i}.root CC.cal ../1769_skim/skim0${i}.dat
  ./r2root /run/media/jing.li/June172018/user/1769CoulTK2/TK${i}.dat ../1769_coulex_root/coulex0${i}.root CC.cal
done
