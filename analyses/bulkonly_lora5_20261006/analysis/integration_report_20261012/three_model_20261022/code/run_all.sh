set -x
L=/mnt/user-data/outputs/s3/logs
python3 comp_s3.py > $L/comp.log 2>&1
python3 myeloid_s3.py sepsis > $L/myeloid_sepsis.log 2>&1
for c in sepsis covid; do python3 c2st_dc.py $c > $L/c2st_$c.log 2>&1; done
for c in sepsis covid; do python3 coexpr.py $c > $L/coexpr_$c.log 2>&1; done
for c in sepsis covid; do python3 lr.py $c > $L/lr_$c.log 2>&1; done
for c in sepsis covid; do python3 claims_rare.py $c > $L/claims_$c.log 2>&1; done
python3 myeloid_s3.py covid > $L/myeloid_covid.log 2>&1
for c in sepsis covid; do python3 integ.py $c > $L/integ_$c.log 2>&1; done
echo ALLDONE
