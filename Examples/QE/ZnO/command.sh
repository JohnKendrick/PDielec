python=$1
shift
params=$*
$python ../../../pdielec $params -matrix ptfe -sigma 5 \
                         -needle 0 0 1 -ellipsoid 0 0 1 0.5 -plate 1 0 0 \
                         -LO 1 0 0  -LO 0 1 0  -LO 0 0 1 \
                         -csv command.csv \
                         -masses program -program qe zno.ph.dynG  $*
