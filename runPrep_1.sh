#!/bin/bash
#Need to modify path, genomes and parameters!

#$ -wd $HOME/lastZ
#$ -o $HOME/lastZ/runPrep.$JOB_ID.out

chainParams="-minScore=1000 -linearGap=medium" #loose for galGal2
lastzParams="C=0 E=30 H=2000 K=3000 L=3000 O=400 M=50" #parameters from UCSC webpage for human-mouse alignments M=50 form cow-pig align in UCSC
TNAME=<sps_name>
QNAME=<sps_name2>

BASEDIR=$(pwd)
LOGDIR=$(pwd)/logs
WRKDIR=$(pwd)/data
TARGET=${WRKDIR}/${TNAME}.2bit
QUERY=${WRKDIR}/${QNAME}.2bit

source /etc/profile.d/modules.sh
module load apps/kent

cd $WRKDIR

ls -ld $TARGET $QUERY

if [ ! -s ${TNAME}.chrom.sizes ]; then
    twoBitInfo ${TARGET} stdout | sort -k2nr > ${TNAME}.chrom.sizes
    rm -fr ${TNAME}PartList ${TNAME}.part.list
    mkdir ${TNAME}PartList
fi
if [ ! -s ${QNAME}.chrom.sizes ]; then
    twoBitInfo ${QUERY} stdout | sort -k2nr > ${QNAME}.chrom.sizes
    rm -fr ${QNAME}PartList ${QNAME}.part.list
    mkdir ${QNAME}PartList
fi

if [ ! -s ${TNAME}.part.list ]; then
    $BASEDIR/scripts/partitionSequence.pl 40010000 10000 ${TARGET} ${TNAME}.chrom.sizes 1 \
        -lstDir ${TNAME}PartList > ${TNAME}.part.list
fi
if [ ! -s ${QNAME}.part.list ]; then
    $BASEDIR/scripts/partitionSequence.pl 20000000 0 ${QUERY} ${QNAME}.chrom.sizes 1 \
        -lstDir ${QNAME}PartList > ${QNAME}.part.list
fi

grep -v PartList ${TNAME}.part.list > target.list
for F in ${TNAME}PartList/*.lst; do
    cat ${F}
done >> target.list

grep -v PartList ${QNAME}.part.list > query.list
for F in ${QNAME}PartList/*.lst; do
    cat ${F}
done >> query.list

$BASEDIR/scripts/constructLiftFile.pl ${TNAME}.chrom.sizes target.list > target.lift
$BASEDIR/scripts/constructLiftFile.pl ${QNAME}.chrom.sizes query.list > query.lift

QS=$(wc -l query.list | cut -f1 -d' ')
TS=$(wc -l target.list | cut -f1 -d' ')
PRODUCT=$(($QS * $TS))

cat <<_EOF_ | sed 's/^##/#/g' > $BASEDIR/scripts/runLastz.sh
#!/bin/bash
##$ -wd $WRKDIR
##$ -o $LOGDIR/runLastz/\$JOB_NAME.\$JOB_ID.\$TASK_ID
source /etc/profile.d/modules.sh
module load apps/lastz apps/kent

B=\$((\$SGE_TASK_ID-1))
T_IDX=\$((\$B / ${QS} + 1))
Q_IDX=\$((\$B % ${QS}))

T=\$(tail -n+\${T_IDX} target.list | head -n1)
Q=\$(tail -n+\${Q_IDX} query.list | head -n1)
FT=\$(basename \$T)
FQ=\$(basename \$Q)

tmpDir=~/localscratch/\${SGE_TASK_ID}.\${FT}
mkdir -p raw psl \${tmpDir}

twoBitToFa \${T} \${tmpDir}/\${FT}.fa
twoBitToFa \${Q} \${tmpDir}/\${FQ}.fa

lastz \${tmpDir}/\${FT}.fa \${tmpDir}/\${FQ}.fa \\
    $lastzParams \\
    > raw/\${FT}.\${FQ}.lav

lavToPsl raw/\${FT}.\${FQ}.lav stdout | \\
    liftUp -type=.psl stdout target.lift error stdin | \\
    liftUp -nohead -pslQ -type=.psl stdout query.lift error stdin | \\
    gzip -c \\
    > psl/\${FT}.\${FQ}.psl.gz

rm -f \${tmpDir}/\${FT}.fa \${tmpDir}/\${FQ}.fa
rmdir --ignore-fail-on-non-empty \${tmpDir}
_EOF_

#HERE!

mkdir -p concat $LOGDIR/runLastz
cat <<EOF | sed 's/^##/#/g' > $BASEDIR/scripts/runConcat.sh
#!/bin/bash
##$ -o $LOGDIR/runConcat.\$JOB_ID.out
##$ -wd $WRKDIR

EOF
perl -pe 's/\:\d+\-\d+//g' target.list > $BASEDIR/data/target_shortNameUnSort
sort -u $BASEDIR/data/target_shortNameUnSort > $BASEDIR/data/target_shortName
for T in `cat target_shortName | sed -e "s#${WRKDIR}/##"`; do
    cat <<EOF >> $BASEDIR/scripts/runConcat.sh
zcat ${WRKDIR}/psl/${T}\:*.psl.gz \\
    | perl -ne '!/^\#/ && print' >> concat/${T}.psl
EOF
done

mkdir -p chain $LOGDIR/runLastz
mkdir -p net $LOGDIR/runLastz
cat <<EOF | sed 's/^##/#/g' > $BASEDIR/scripts/runChainNet.sh
#!/bin/bash
##$ -o $LOGDIR/runChain.\$JOB_ID.\$TASK_ID.out
##$ -wd $WRKDIR
#$ -l h_vmem=96G

source /etc/profile.d/modules.sh
module load apps/kent
EOF

L=$(wc -l ${WRKDIR}/${TNAME}.chrom.sizes | cut -f1 -d' ')

#modify it and use FIND instead
#find ./chain -name "*.chain" | \\
cat <<EOF >> $BASEDIR/scripts/runChainNet.sh

find ${WRKDIR}/concat -name "*.psl" >${WRKDIR}/find_list 

B=\$((\$SGE_TASK_ID-1))
IDX=\$((\$B % ${L}))

T=\$(tail -n+\${IDX} find_list | head -n1)
FT=\$(basename \$T)

DIR=${WRKDIR}/concat
FILE=${DIR}/${FT}
  
axtChain -psl -verbose=0 ${chainParams} \${FILE} ${TARGET} ${QUERY} stdout \
| chainAntiRepeat ${TARGET} ${QUERY} stdin \${FILE}.chain 
chainSort \${FILE}.chain stdout | chainPreNet stdin ${WRKDIR}/${TNAME}.chrom.sizes ${WRKDIR}/${QNAME}.chrom.sizes stdout \
| chainNet stdin ${WRKDIR}/${TNAME}.chrom.sizes ${WRKDIR}/${QNAME}.chrom.sizes stdout /dev/null \
| netSyntenic stdin \${FILE}.net

EOF
cat <<EOF >>$BASEDIR/scripts/runChainNet.sh
mv ${WRKDIR}/concat/*.chain ${WRKDIR}/chain/
mv ${WRKDIR}/concat/*.net ${WRKDIR}/net/
EOF

cat <<EOF
Parameters:
chains/nets: ${chainParams}
lastZ: ${lastzParams}
Preparation is complete.  Run:

  qsub -t 1-${PRODUCT} $BASEDIR/scripts/runLastz.sh

After the alignments are completed, run:

  qsub $BASEDIR/scripts/runConcat.sh

After Concat is done, run:

  qsub -t 1-${L} $BASEDIR/scripts/runChainNet.sh

EOF


