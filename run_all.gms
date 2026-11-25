* Runs the tri-stage model for the three comfort levels used in Tables 10-11 and Figures 6-12.
* Open this file in GAMS Studio (or: gams run_all.gms) - it writes results_0.01/0.5/0.99.csv and .gdx
* into the working folder; then run make_figures.py.
$call gams TriStage_CloudFogEdge.gms --ALPHA=0.01 lo=2
$if errorlevel 1 $abort 'run for ALPHA=0.01 failed - see TriStage_CloudFogEdge.lst'
$call gams TriStage_CloudFogEdge.gms --ALPHA=0.5 lo=2
$if errorlevel 1 $abort 'run for ALPHA=0.5 failed - see TriStage_CloudFogEdge.lst'
$call gams TriStage_CloudFogEdge.gms --ALPHA=0.99 lo=2
$if errorlevel 1 $abort 'run for ALPHA=0.99 failed - see TriStage_CloudFogEdge.lst'
$log results_0.01, results_0.5 and results_0.99 (.csv and .gdx) written
