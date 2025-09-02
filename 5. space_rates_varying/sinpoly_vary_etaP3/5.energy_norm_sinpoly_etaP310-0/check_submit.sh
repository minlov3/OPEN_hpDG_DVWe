#!/bin/bash
#SBATCH -J hp_sinsine_t_timeh
#SBATCH -p node
#SBATCH -N 1
#SBATCH -n 40
##以上第1-5行为作业提交相关参数，具体可参考《Slurm作业调度系统简介》

MATLAB_HOME=/share/apps/Polyspace/R2020a
export PATH=$MATLAB_HOME/bin:$PATH
##以上第8-9行为MATLAB运行相关环境变量，保持不变

date
matlab -nodesktop -nosplash -nodisplay -r Check_with_energy_viscous_Wave_PDE_2D_sinpoly>>  space_h_eta.log
date
##第13行为调用MATLAB对执行源代码的命令，其中“-r”后为MATLAB源码文件名，本例中为“ceshi”，注意不带“.m”后缀
