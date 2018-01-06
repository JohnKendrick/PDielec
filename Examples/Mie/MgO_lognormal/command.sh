python=$1
shift
params=$*
$python ../../../pdielec $params -program castep phonon -vmin 300 -vmax 800 -sphere -vf 0.1 -method mie  -size 1.0 0.1 -sigma 5 -csv command.csv $*
