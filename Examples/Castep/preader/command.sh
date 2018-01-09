python=$1
shift
params=$*
$python ../../../preader $params -program castep -masses isotopic ../AsparticAcid/phonon.castep  ../MgO/phonon.castep  ../Na2SO42/phonon.castep 
