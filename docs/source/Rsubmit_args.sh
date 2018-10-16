#!/bin/bash -l
#PBS -l nodes=1:ppn=3,walltime=03:00:00

R --no-save --args ${I} < ${script}.R > ${script}.Rout
