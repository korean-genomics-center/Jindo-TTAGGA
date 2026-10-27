#!/bin/bash
set -uo pipefail
set +u
conda activate TTT
set -u
python -c "import h5py" 2>/dev/null || { echo "FATAL: h5py 없음"; exit 1; }
FA=$(sed -n "${SGE_TASK_ID}p" ${JINDO_ROOT}/Analysis/Population_20Jindo/Reanalysis_Combined/ROS_RepeatMask/task.list)
BASE=$(basename $FA .fa)
echo "=== ROS RepeatMask $BASE : $(date) on $(hostname) ==="
OUTF=${JINDO_ROOT}/Analysis/Population_20Jindo/Reanalysis_Combined/ROS_RepeatMask/out/${BASE}.fa.out
[ -s "$OUTF" ] && { echo "SKIP $BASE"; exit 0; }
${JINDO_ROOT}/Resources/Tools/RepeatAnnotation_Masker_Modeler/RepeatMasker/RepeatMasker_v4.2.3/RepeatMasker/RepeatMasker -e rmblast -pa 2 -a -s -species dog -xsmall -dir ${JINDO_ROOT}/Analysis/Population_20Jindo/Reanalysis_Combined/ROS_RepeatMask/out $FA
rc=$?
if [ $rc -ne 0 ] || [ ! -s "$OUTF" ]; then echo "FATAL: $BASE (rc=$rc)"; exit 1; fi
echo "=== DONE $BASE : $(wc -l < $OUTF) lines : $(date) ==="
