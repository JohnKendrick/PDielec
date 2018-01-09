python=$1
shift
params=$*
$python ../../../pdielec $params -matrix ptfe -sigma 5 \
                         -needle 0 0 1 -ellipsoid 0 0 1 0.5 -plate 1 0 0 \
                         -LO 1 1 1 \
                         -csv command1.csv \
                         -masses program -program qe phonon.dynG  $*
$python ../../../pdielec $params -matrix ptfe -sigma 5 \
                         -sphere -masses isotope -print \
                         -csv command2.csv \
                         -program qe phonon.dynG  $*
$python ../../../pdielec $params -matrix ptfe -sigma 5 \
                         -sphere -masses isotope -print \
                         -mass H 2.01410178 \
                         -csv command3.csv \
                         -program qe phonon.dynG  $*
cat command1.csv command2.csv command3.csv > command.csv
rm -f command1.csv command2.csv command3.csv
