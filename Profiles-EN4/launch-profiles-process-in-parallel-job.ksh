#!/bin/sh
#SBATCH --job-name=profiles-process 
#SBATCH --ntasks=145
#SBATCH --ntasks-per-node=40       
#SBATCH --hint=nomultithread       
#SBATCH --time=00:10:00            
#SBATCH --output=profiles-process%j.out  
#SBATCH --error=profiles-process%j.out  
#SBATCH --account=egi@cpu

NB_NPROC=145 ##(4 var * 475 dates)

runcode() { srun --mpi=pmi2 -m cyclic -n $@ ; }
liste=''

for k in $(seq 1 $NB_NPROC); do
    kk=`expr $k + 1`
    cp tmp_script_profile.ksh script_profile_${k}.ksh
    chmod +x script_profile_${k}.ksh
    sed -i "s/PROF/$k/g" script_profile_${k}.ksh
    liste="$liste /gpfswork/rech/egi/rote001/git/diags-CMEMS-on-jean-zay/Profiles-EN4/script_profile_${k}.ksh"

done

runcode $NB_NPROC /gpfswork/rech/egi/rote001/git/DMONTOOLS/MPI_TOOLS/mpi_shell $liste
