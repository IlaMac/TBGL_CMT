#!/bin/bash


############### SET MAIN FOLDERS ###############
BASEDIR=${HOME}/TBGL_CMT
SCRIPT_DIR=${BASEDIR}/launch_scripts

cd /tmp/

if [ ! -d ./SOutput_x_ilaria ]; then
   mkdir -p Output_x_ilaria
fi
#############################################

time_limit="7-00:00:00"

LLIST="64"


RESTART=0
#RESTART=0-> Restart from scratch
#RESTART=1-> Restart from interrupted run
#RESTART=2-> Restart from the previous final scenario

############ Parameters of the Hamiltonian ---> HP_init.txt in a directory whose name contains the main parameters values##################
H_K=-1
H_lambda=0.1
H_e=0
H_h=1
H_blow=3.0
H_bhigh=8.0

ExtField=1	# 0: no Field; 1: Field
if [[ "${ExtField}" == "1" ]]; then
	H_fx=0.
	H_fy=0.015625 #1/64
fi

if [[ "${ExtField}" == "0" ]]; then
	H_fx=0
	H_fy=0
fi

H_init=4
#If H_init=0: ground state K>0: \th_12=+pi/2; \gamma=pi/2
#H_init=1: ground state K<0; lambda>0; \th_12=0; gamma=C_PI/3;
#H_init=2: ground state K<0; lambda>0; gamma=C_PI
#H_init=3: London anisotropic
#H_init=4: random

H_london=0
#H_london=0 density update with total density constraint
#H_london=1 London approximation

############ Parameters for the Monte Carlo simulations --> MC_init.txt#####################

Nmisu=300000
ntau=32
ntransient=100000
freq_autosave=1000
l_box=1
theta_box=3.141592653
A_box=0.1
#meas_corr=0 #compute ot not the correlation functions

#############################################

for L in $LLIST; do

############Creation of the output folder and of the two files of initialization####################

cd ${BASEDIR}/Output_TBGL_CMT

if [[ "${ExtField}" == "1" ]]; then
	if [ ! -d ./Sfx_${H_fx}_fy_${H_fy} ]; then
   		mkdir -p fx_${H_fx}_fy_${H_fy}
	fi
	cd fx_${H_fx}_fy_${H_fy}
fi

if [ ! -d ./Slambda_${H_lambda} ]; then
   mkdir -p lambda_${H_lambda}
fi

cd lambda_${H_lambda}

if [ ! -d ./SK_${H_K} ]; then
   mkdir -p K_${H_K}
fi

cd K_${H_K}

if [ ! -d ./Se_${H_e} ]; then
   mkdir -p e_${H_e}
fi

cd e_${H_e}

if [ ! -d ./Sh_${H_h} ]; then
   mkdir -p h_${H_h}
fi

cd h_${H_h}

if [ ! -d ./SL${L}_lambda${H_lambda}_K${H_K}_e${H_e}_h${H_h}_bmin${H_blow}_bmax${H_bhigh}__init${H_init}_london${H_london} ]; then
   mkdir -p L${L}_lambda${H_lambda}_K${H_K}_e${H_e}_h${H_h}_bmin${H_blow}_bmax${H_bhigh}_init${H_init}_london${H_london}
fi

if [[ "${ExtField}" == "1" ]]; then
  OUTPUT=${BASEDIR}/Output_TBGL_CMT/fx_${H_fx}_fy_${H_fy}/lambda_${H_lambda}/K_${H_K}/e_${H_e}/h_${H_h}/L${L}_lambda${H_lambda}_K${H_K}_e${H_e}_h${H_h}_bmin${H_blow}_bmax${H_bhigh}_init${H_init}_london${H_london}
fi

if [[ "${ExtField}" == "0" ]]; then
  OUTPUT=${BASEDIR}/Output_TBGL_CMT/lambda_${H_lambda}/K_${H_K}/e_${H_e}/h_${H_h}/L${L}_lambda${H_lambda}_K${H_K}_e${H_e}_h${H_h}_bmin${H_blow}_bmax${H_bhigh}_init${H_init}_london${H_london}
fi

#############################################
cd /tmp/Output_x_ilaria

if [[ "${ExtField}" == "1" ]]; then
	if [ ! -d ./Sfx_${H_fx}_fy_${H_fy} ]; then
   		mkdir -p fx_${H_fx}_fy_${H_fy}
	fi
	cd fx_${H_fx}_fy_${H_fy}
fi


if [ ! -d ./Slambda_${H_lambda} ]; then
   mkdir -p lambda_${H_lambda}
fi

cd lambda_${H_lambda}

if [ ! -d ./SK_${H_K} ]; then
   mkdir -p K_${H_K}
fi

cd K_${H_K}

if [ ! -d ./Se_${H_e} ]; then
   mkdir -p e_${H_e}
fi

cd e_${H_e}

if [ ! -d ./Sh_${H_h} ]; then
   mkdir -p h_${H_h}
fi

cd h_${H_h}

if [ ! -d ./SL${L}_lambda${H_lambda}_K${H_K}_e${H_e}_h${H_h}_bmin${H_blow}_bmax${H_bhigh}_init${H_init}_london${H_london} ]; then
   mkdir -p L${L}_lambda${H_lambda}_K${H_K}_e${H_e}_h${H_h}_bmin${H_blow}_bmax${H_bhigh}_init${H_init}_london${H_london}
fi

if [[ "${ExtField}" == "1" ]]; then
  OUTPUT_TEMP=/tmp/Output_x_ilaria/fx_${H_fx}_fy_${H_fy}/lambda_${H_lambda}/K_${H_K}/e_${H_e}/h_${H_h}/L${L}_lambda${H_lambda}_K${H_K}_e${H_e}_h${H_h}_bmin${H_blow}_bmax${H_bhigh}_init${H_init}_london${H_london}
fi

if [[ "${ExtField}" == "0" ]]; then
  OUTPUT_TEMP=/tmp/Output_x_ilaria/lambda_${H_lambda}/K_${H_K}/e_${H_e}/h_${H_h}/L${L}_lambda${H_lambda}_K${H_K}_e${H_e}_h${H_h}_bmin${H_blow}_bmax${H_bhigh}_init${H_init}_london${H_london}
fi
#############################################

cd ${OUTPUT}

#THE ORDER OF WRITING DOES MATTER
echo $H_K > HP_init.txt
echo $H_lambda >> HP_init.txt
echo $H_e >> HP_init.txt
echo $H_h >> HP_init.txt
echo $H_blow >> HP_init.txt
echo $H_bhigh >> HP_init.txt
echo $H_fx >> HP_init.txt
echo $H_fy >> HP_init.txt
echo $H_init >> HP_init.txt
echo $H_london>> HP_init.txt

#THE ORDER OF WRITING DOES MATTER
echo $Nmisu > MC_init.txt
echo $ntau >> MC_init.txt
echo $ntransient >> MC_init.txt
echo $freq_autosave >> MC_init.txt
echo $l_box >> MC_init.txt
echo $theta_box >> MC_init.txt
echo $A_box >> MC_init.txt
#echo $meas_corr >> MC_init.txt

#################Creation of the submit_runs script#########################

jobname="L${L}_lambda${H_lambda}_K${H_K}_e${H_e}_h${H_h}_bmin${H_blow}_bmax${H_bhigh}_init${H_init}_london${H_london}_fx${H_fx}_fy${H_fy}"
nnodes=2
ntasks=64 #parallel tempering over ntasks temperatures

#I create ntasks folder: one for each rank.

cd ${OUTPUT}

for ((rank=0; rank<${ntasks}; rank++)); do

if [ ! -d ./Sbeta_${rank} ]; then
   mkdir -p beta_${rank}
fi

done

cd ${OUTPUT_TEMP}

for ((rank=0; rank<${ntasks}; rank++)); do

if [ ! -d ./Sbeta_${rank} ]; then
   mkdir -p beta_${rank}
fi

done

cd ${SCRIPT_DIR}
DIR_PAR="${OUTPUT}"
DIR_PAR_TEMP="${OUTPUT_TEMP}"

#SEED= If I want to repeat exactly a simulation I could initialize the random number generator exactly at the same way

EXECUTE_DIR="../build/release-conan"

#SBATCH --nodes=${nnodes}               # Number of nodes

echo "#!/bin/bash
#SBATCH --job-name=${jobname}          # Name of the job
#SBATCH --time=${time_limit}               # Allocation time
#SBATCH --mem-per-cpu=2000              # Memory per allocated cpu
#SBATCH --nodes=${nnodes}               # Number of nodes
#SBATCH --ntasks=${ntasks}
#SBATCH --output=${DIR_PAR}/logs/log_${jobname}.o
#SBATCH --error=${DIR_PAR}/logs/log_${jobname}.e

srun ${EXECUTE_DIR}/CMT -L ${L} --outdir ${DIR_PAR} --tempdir ${DIR_PAR} -r ${RESTART} -s -1 &> ${DIR_PAR}/logs/log_${jobname}.o

" >  submit_run

#Submission of the work --> sbatch submit_runs

mkdir -p ${DIR_PAR}/logs

sbatch submit_run

done
