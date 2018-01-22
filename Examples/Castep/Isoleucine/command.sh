python=$1
shift
params=$*
$python ../../../phonana -program castep phonon.castep -csv command.csv  $*
