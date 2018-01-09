python=$1
shift
params=$*
$python ../../../preader $params -program phonopy vasp -eckart ../ZnO/OUTCAR ../Na2SO42/OUTCAR

