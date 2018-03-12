python=$1
shift
params=$*
$python ../../../pdielec $params -eckart -neutral -masses program -mf 0.1 -matrix nujol -vmax 350 -sigma 5 -sphere -plate 1 0 0 -plate 0 1 0  -plate 0 0 1  -method maxwell phonon -csv command.csv  $*
