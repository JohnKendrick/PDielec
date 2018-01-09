python=$1
shift
params=$*
$python ../../../pdielec -size 0 -size 1 -size 3 $params -method ap -method bruggeman -matrix ptfe -sigma 5 \
                 -ellipsoid 0 0 1 0.5 \
                 -csv command1.csv \
                 -masses program -program abinit BaTiO3.out $* 

$python ../../../pdielec -size 0 -size 1 -size 3 $params -method ap -method bruggeman -matrix ptfe -sigma 5 \
                 -plate 1 0 0 \
                 -csv command2.csv \
                 -masses average \
                 -print \
                 -program abinit BaTiO3.out $* 
cat command1.csv command2.csv > command.csv
rm -f command1.csv command2.csv
