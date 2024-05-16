## This document contains notes for polishing the t2t primates assemblies

### Check read coverages for all data

Location of data in this spreadsheet: https://docs.google.com/spreadsheets/d/17HVGUJ7cvf8SvxkLNVg4cGcIvDTaDCJvnGf56CT27d4/edit#gid=437306236

```
/private/groups/patenlab/mira/t2t_primates_polishing/read_coverages
```

```sh
#SBATCH --job-name=single-node-run
#SBATCH --cpus-per-task=8
#SBATCH --threads-per-core=1
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=16gb
#SBATCH --time=7-00:00
#SBATCH --partition=long
#SBATCH --output=slurm_logs/submission_%x_%j_%A_%a.log


SAMPLE_ID=$(awk -F ',' -v task_id=${SLURM_ARRAY_TASK_ID} 'NR>1 && NR==task_id+1 {print $1}' "${SAMPLE_CSV}")

# Ensure a sample ID is obtained
if [ -z "${SAMPLE_ID}" ]; then
    echo "Error: Failed to retrieve a valid sample ID. Exiting."
    exit 1
fi
```
