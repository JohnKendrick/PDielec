python=$1
shift
params=$*
$python ../../../pdielec $params -method ap -method maxwell -matrix ptfe -sigma 5 \
                 -plate 0 0 1 -ellipsoid 0 0 1 0.5 -plate 1 0 0 \
                 -csv command.csv \
                 -program abinit AlAs.out $*
