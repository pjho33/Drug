#!/bin/bash
# Run control complex MD simulation with OpenMM

cd /home/pjho3/projects/Drug/final_complex/controlcomplex/openmm

echo "=========================================="
echo "Starting Control Complex MD Simulation"
echo "=========================================="
echo ""

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate Drug-MD

# Set variables
init="step5_input"
equi_prefix="step6.%d_equilibration"
prod_prefix="step7_production"
prod_step="step7"

echo "Step 1: Running Equilibration (6 steps)"
echo "------------------------------------------"

# Equilibration steps 1-6
for cnt in {1..6}; do
    pcnt=$((cnt - 1))
    istep=$(printf "${equi_prefix}" ${cnt})
    pstep=$(printf "${equi_prefix}" ${pcnt})
    
    echo ""
    echo "Running ${istep}..."
    
    if [ ${cnt} -eq 1 ]; then
        # First equilibration step
        python -u openmm_run.py \
            -i ${istep}.inp \
            -t toppar.str \
            -p ${init}.psf \
            -c ${init}.crd \
            -b sysinfo.dat \
            -orst ${istep}.rst \
            -odcd ${istep}.dcd \
            > ${istep}.out 2>&1
    else
        # Subsequent equilibration steps
        python -u openmm_run.py \
            -i ${istep}.inp \
            -t toppar.str \
            -p ${init}.psf \
            -c ${init}.crd \
            -irst ${pstep}.rst \
            -orst ${istep}.rst \
            -odcd ${istep}.dcd \
            > ${istep}.out 2>&1
    fi
    
    if [ $? -eq 0 ]; then
        echo "✅ ${istep} completed"
    else
        echo "❌ ${istep} failed"
        echo "Check ${istep}.out for errors"
        exit 1
    fi
done

echo ""
echo "=========================================="
echo "✅ Equilibration completed!"
echo "=========================================="
echo ""
echo "Step 2: Running Production (100ns)"
echo "------------------------------------------"

# Production run
cnt=1
pstep=$(printf "${equi_prefix}" 6)

echo ""
echo "Running ${prod_step}_${cnt}..."

python -u openmm_run.py \
    -i ${prod_prefix}.inp \
    -t toppar.str \
    -p ${init}.psf \
    -c ${init}.crd \
    -irst ${pstep}.rst \
    -orst ${prod_step}_${cnt}.rst \
    -odcd ${prod_step}_${cnt}.dcd \
    > ${prod_step}_${cnt}.out 2>&1

if [ $? -eq 0 ]; then
    echo "✅ ${prod_step}_${cnt} completed"
else
    echo "❌ ${prod_step}_${cnt} failed"
    echo "Check ${prod_step}_${cnt}.out for errors"
    exit 1
fi

echo ""
echo "=========================================="
echo "✅ Production run completed!"
echo "=========================================="
echo ""
echo "Output files:"
echo "  Equilibration: step6.*.dcd, step6.*.rst"
echo "  Production: step7_1.dcd, step7_1.rst"
echo ""
