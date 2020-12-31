#!/bin/bash
# compile_run.sh [ TIMES ]

if [ $# -lt 2 ]; then
  TIMES=1
else
  TIMES="$2"
fi

mpic++ flpenum.cpp -o flpenum.out -w

for i in $(seq $TIMES)
do 
  mpirun -np $1 --oversubscribe ./flpenum.out $3 | tee -a out.log
  echo "-------------------" | tee -a out.log
done

echo "===================" >> out.log