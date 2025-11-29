#!/bin/bash
# Monitor domain decomposition test jobs

echo "========================================"
echo "Domain Decomposition Test Jobs Status"
echo "========================================"
echo "Time: $(date '+%Y-%m-%d %H:%M:%S')"
echo ""

cd /scratch/gpfs/hshao/C_Bible/FFT_Hermitian/fftw_openmp/mpi_mock/hermitian_3d_matrix_production/tests/domain_decomposition/slurm_jobs

if [ ! -f test_jobs.txt ]; then
    echo "No jobs found. Run submit_all_tests.sh first."
    exit 1
fi

printf "%-12s %-10s %-8s %-8s %-8s %-15s %-10s\n" "Job ID" "Test" "N" "Ranks" "Nodes" "Status" "Time"
printf "%-12s %-10s %-8s %-8s %-8s %-15s %-10s\n" "------" "----" "-" "-----" "-----" "------" "----"

while IFS=' ' read -r job_id test_name N ranks nodes; do
    if [[ $job_id == \#* ]] || [[ -z $job_id ]]; then
        continue
    fi
    
    # Get job status
    status=$(squeue -j $job_id -h -o "%T" 2>/dev/null)
    
    if [ -z "$status" ]; then
        # Job not in queue, check if completed
        out_file="slurm-${test_name}.out"
        if [ -f "$out_file" ] && grep -q "Test complete" "$out_file"; then
            status="COMPLETED"
            # Try to extract time from log
            time=$(grep -i "elapsed\|time" "$out_file" | tail -1 | grep -oP '[0-9]+\.[0-9]+' | head -1 || echo "N/A")
        else
            status="FAILED/UNKNOWN"
            time="N/A"
        fi
    else
        time=$(squeue -j $job_id -h -o "%M" 2>/dev/null)
    fi
    
    printf "%-12s %-10s %-8s %-8s %-8s %-15s %-10s\n" "$job_id" "$test_name" "$N" "$ranks" "$nodes" "$status" "$time"
done < test_jobs.txt

echo ""
echo "To check detailed output:"
echo "  cat slurm-*.out"
echo ""
echo "To check test results:"
echo "  cd ../results && ls -lh"
echo ""

