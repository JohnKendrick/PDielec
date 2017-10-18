python=$1
shift
params=$*
$python ../../../preader $params -program qe ../Cocaine/phonon.dynG  ../Na2SO42/Na2SO42.dynG ../ZnO/zno.ph.dynG
