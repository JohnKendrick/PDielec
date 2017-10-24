python=$1
shift
params=$*
$python ../../../pdielec $params  -size 0 -size 0.1 -size 1 -size 3 -vmin 300 -vmax 800 -sphere -method maxwell   -sigma 10 -csv command1.csv phonon  $*
$python ../../../pdielec $params  -size 0 -size 0.1 -size 1 -size 3 -vmin 300 -vmax 800 -sphere -method bruggeman -sigma 10 -csv command2.csv phonon  $*
cat command1.csv command2.csv > command.csv
rm command1.csv command2.csv
