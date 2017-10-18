python=$1
shift
params=$*
$python ../../../pdielec $params -matrix ptfe -vmax 4000 -sigma 2 -sphere -plate 1 0 0 -plate 0 0 1  -ellipsoid 0 0 1  0.5 -method ap -method maxwell phonon -csv command.csv  $*
