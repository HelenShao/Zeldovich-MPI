#!/bin/bash
# Submit all domain decomposition test jobs via SLURM

echo "========================================"
echo "Submitting Domain Decomposition Test Jobs"
echo "========================================"
echo "Date: $(date '+%Y-%m-%d %H:%M:%S')"
echo ""

cd /scratch/gpfs/hshao/C_Bible/FFT_Hermitian/fftw_openmp/mpi_mock/hermitian_3d_matrix_production/tests/domain_decomposition/slurm_jobs

# Create job tracking file
echo "# Job_ID Test_Name N Ranks Nodes" > test_jobs.txt

# Function to calculate nodes needed (assuming 32 ranks per node)
calculate_nodes() {
    local ranks=$1
    local ranks_per_node=32
    local nodes=$(( (ranks + ranks_per_node - 1) / ranks_per_node ))
    echo $nodes
}

# Function to create and submit a test job
submit_test() {
    local test_name=$1
    local N=$2
    local ranks=$3
    
    local nodes=$(calculate_nodes $ranks)
    
    # Create SLURM script from template
    local slurm_file="${test_name}.slurm"
    
    # Replace template variables
    sed -e "s|\${TEST_NAME}|${test_name}|g" \
        -e "s|\${N}|${N}|g" \
        -e "s|\${RANKS}|${ranks}|g" \
        -e "s|\${NODES}|${nodes}|g" \
        template.slurm > "$slurm_file"
    
    # Submit job
    JOB_ID=$(sbatch "$slurm_file" | awk '{print $4}')
    
    if [ ! -z "$JOB_ID" ]; then
        echo "Submitted: $test_name (N=$N, Ranks=$ranks, Nodes=$nodes, Job ID: $JOB_ID)"
        echo "$JOB_ID $test_name $N $ranks $nodes" >> test_jobs.txt
    else
        echo "Failed to submit: $test_name"
    fi
}

# Exact Division Tests (from TEST_PLAN_DOMAIN_DECOMPOSITION.md)
echo "Submitting Exact Division Tests..."
submit_test "E1" 1024 16
submit_test "E2" 512 8
submit_test "E3" 256 4
submit_test "E6" 1000 16

# Remainder Tests (from TEST_PLAN_DOMAIN_DECOMPOSITION.md)
echo ""
echo "Submitting Remainder Tests..."
submit_test "R1" 1001 16
submit_test "R2" 1025 16
submit_test "R3" 1003 16
submit_test "R4" 1024 13
submit_test "R5" 1000 13
submit_test "R6" 1001 13

# Non-Power-of-2 Tests (from TEST_PLAN_DOMAIN_DECOMPOSITION.md)
echo ""
echo "Submitting Non-Power-of-2 Tests..."
submit_test "NP1" 175 16
submit_test "NP2" 175 25
submit_test "NP3" 175 7
submit_test "NP4" 175 35
submit_test "NP5" 175 8
submit_test "NP6" 175 13

# Note: These tests use the exact ranks specified in TEST_PLAN_DOMAIN_DECOMPOSITION.md
# Some tests may use fewer ranks than the minimum required (N/2 + 1), which will cause
# the code to fail with "Too few ranks" error. This is expected for testing purposes.

echo ""
echo "========================================"
echo "All jobs submitted!"
echo "Job IDs saved to: test_jobs.txt"
echo ""
echo "To monitor jobs:"
echo "  ./monitor_test_jobs.sh"
echo ""
echo "To collect results (after completion):"
echo "  cd .. && ./run_all_tests.sh"
echo ""

