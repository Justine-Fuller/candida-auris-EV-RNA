# helper (not a job): known-answer test job, then array, then summary job with dependencies
B=/scratch/flljus002/Phase3D_forensic_audit/v2_read_census_2026-10-09
t=$(sbatch --parsable -p ada -c 1 --mem=2G -t 00:10:00 -J census_test -o $B/logs/census_test_%j.out --wrap "/opt/exp_soft/miniconda3-py3.12/bin/python3 $B/scripts/test_census.py $B/scripts $B/run/test | tee $B/run/test_result.txt; grep -q ALL_CENSUS_TESTS_PASS $B/run/test_result.txt")
a=$(sbatch --parsable --dependency=afterok:$t $B/scripts/census.sbatch)
s=$(sbatch --parsable --dependency=afterok:$a -p ada -c 1 --mem=8G -t 00:30:00 -J census_summary -o $B/logs/census_summary_%j.out --wrap "/opt/exp_soft/miniconda3-py3.12/bin/python3 $B/scripts/summarize_census.py $B/run/counts $B/manifest.tsv $B/run/summary | tee $B/run/summary_log.txt; cd $B/run && md5sum counts/* summary/* > outputs.md5 && chmod -R a-w counts summary")
echo "test=$t array=$a summary=$s"
