python=$1
shift
params=$*
$python ../../../pdielec $params -program castep phonon -vmin 300 -vmax 800 -sphere -vf 0.01 -vf 0.1 -vf 0.2 -method maxwell -method mie  -size 0 -size 0.1 -size 1.2 -sigma 5 -csv command.csv $*
