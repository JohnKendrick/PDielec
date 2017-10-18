python=$1
shift
params=$*
$python ../../../preader $params -program vasp -eckart ../ZnO/OUTCAR ../Na2SO42/OUTCAR ../F-Apatite/OUTCAR

