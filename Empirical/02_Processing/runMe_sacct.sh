#!/bin/env bash

ROOTDIR="/scratch/gas0042/nanopore_benchmark/real_data"

SLURM_DIR="${ROOTDIR}/scripts/logs/slurm_jobid"

JOBID_FILES=(
        "${SLURM_DIR}/flye_slurm_jobid.csv"
        "${SLURM_DIR}/spades_short_slurm_jobid.csv"
        "${SLURM_DIR}/abyss_short_slurm_jobid.csv"
        "${SLURM_DIR}/spades_hybrid_slurm_jobid.csv"
        "${SLURM_DIR}/polypolish_slurm_jobid.csv"
        # "${SLURM_DIR}/abyss-hybrid_slurm_jobid.csv"
        # "${SLURM_DIR}/hifiasm_slurm_jobid.csv"
        # "${SLURM_DIR}/canu_slurm_jobid.csv"
        # "${SLURM_DIR}/pilon_slurm_jobid.csv"
        # "${SLURM_DIR}/masurca_slurm_jobid.csv"
        # "${SLURM_DIR}/raven_slurm_jobid.csv"
        # "${SLURM_DIR}/racon_slurm_jobid.csv"
        # "${SLURM_DIR}/racon-mp2_slurm_jobid.csv"
        # "${SLURM_DIR}/nextdenovo_slurm_jobid.csv"
        # "${SLURM_DIR}/miniasm_slurm_jobid.csv"
        # "${SLURM_DIR}/shasta_slurm_jobid.csv"
)

for JOBID_FILE in "${JOBID_FILES[@]}"; do
    
    # Check if the job ID file exists and is not empty
    if [[ ! -s $JOBID_FILE ]]; then
        echo "Job ID file $JOBID_FILE does not exist or is empty. Skipping."
        continue
        # exit 1
    fi

    JOBID_LIST=$(awk -F','  -v 'ORS=,' '{print $2}' $JOBID_FILE | sed -e "s/JobID,//")

    sacct --format=JobID,JobName%100,Start,End,Elapsed,ReqMem,MaxRSS,ReqCPUS,NCPUS,CPUTime,CPUTimeRAW,Partition,ExitCode,State -P --delimiter="," -j $JOBID_LIST > ${JOBID_FILE%_jobid.csv}_sacct.csv
done

echo "JobID,JobName,Start,End,Elapsed,ReqMem,MaxRSS,ReqCPUS,NCPUS,CPUTime,CPUTimeRAW,Partition,ExitCode,State" > ${SLURM_DIR}/all_sacct_real.csv

cat ${SLURM_DIR}/*slurm_sacct.csv | sed "/JobID,JobName,Start,End,Elapsed,ReqMem,MaxRSS,ReqCPUS,NCPUS,CPUTime,CPUTimeRAW,Partition,ExitCode,State/d" >> ${SLURM_DIR}/all_sacct_real.csv