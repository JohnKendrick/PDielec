python=$1
shift
params=$*
$python ../../../pdielec $params  -size 0 -size 1 -size 2 -matrix ptfe -vmin 0 -vmax 1000 -mf 0.01 -method ap -method maxwell -sphere -needle 0 0 1 -csv command.csv  -sigma 5 .   $*
