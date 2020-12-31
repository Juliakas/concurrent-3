#!/bin/sh
#SBATCH -p short
#SBATCH -n16
#SBATCH -C beta

mpic++ flpenum.cpp -o flpenum.out -w

for i in $(seq $1)
do 
  mpirun -np 16 --oversubscribe ./flpenum.out $2 | tee -a out.log
  echo "-------------------" | tee -a out.log
done

echo "===================" | tee -a out.log