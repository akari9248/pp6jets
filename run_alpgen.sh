#!/bin/bash

# Initialize CMS environment
home="/afs/cern.ch/user/s/shuangyu/"
cd ${home}"CMSSW_10_6_28_patch1"
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
cd -

# Check if chunk number parameter is provided
if [ -z "$1" ]; then
    echo "Error: Please provide chunk number as parameter"
    echo "Usage: $0 <chunk_number>"
    exit 1
fi

CHUNK_NUM=$1
WORK_DIR="alpgen/Njetwork"
OUTPUT_NAME="chunk${CHUNK_NUM}"

# Generate configuration file for step 1
ISEED1=$((RANDOM%90000+10000))
ISEED2=$((RANDOM%90000+10000))

cat > config_step1.txt << EOF
1
"${OUTPUT_NAME}"
0
20000, 1000
20000000
ih2 1
ebeam 6500
ndns 315200
iqopt 1
njets 6
ptjmin 15
etajmax 5
drjmin 0.3
ilhe 1
iseed1 $ISEED1
iseed2 $ISEED2
EOF

ls 
# Generate configuration file for step 2
ISEED3=$((RANDOM%90000+10000))
ISEED4=$((RANDOM%90000+10000))

cat > config_step2.txt << EOF
2
"${OUTPUT_NAME}"
iseed3 $ISEED3
iseed4 $ISEED4
EOF

echo "Configuration files generated with random seeds:"
echo "Step 1: $ISEED1, $ISEED2"
echo "Step 2: $ISEED3, $ISEED4"

# Check if working directory exists
if [ ! -d "$WORK_DIR" ]; then
    echo "Error: Working directory '$WORK_DIR' does not exist"
    exit 1
fi

# Change to working directory
cd "$WORK_DIR" || { echo "Cannot change to directory $WORK_DIR"; exit 1; }

# Clean up previous files
echo "Cleaning up previous files..."
rm -f "${OUTPUT_NAME}"* config_step*.txt

# Move configuration files to current directory
mv ../../config_step*.txt .

# Clean and compile
echo "Compiling generator..."
# make cleanall
# make gen

# Check if compilation was successful
if [ ! -f "Njetgen" ]; then
    echo "Error: Njetgen compilation failed"
    exit 1
fi

# Execute step 1
echo "Running step 1..."
./Njetgen < config_step1.txt

# Execute step 2
echo "Running step 2..."
./Njetgen < config_step2.txt

# Check if final output exists
if [ ! -f "${OUTPUT_NAME}.lhe" ]; then
    echo "Error: ${OUTPUT_NAME}.lhe file not generated"
    exit 1
fi

# Move LHE file to parent directory
mv "${OUTPUT_NAME}.lhe" ../..

echo "Completed! LHE file generated: ../../${OUTPUT_NAME}.lhe"
echo "Temporary files cleaned up"

# Return to original directory
cd - > /dev/null
ls
