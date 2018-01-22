python=$1
shift
params=$*
$python ../../../phonana $params -csv command.csv -radius Ba 0.1 -program abinit BaTiO3.out $* 
