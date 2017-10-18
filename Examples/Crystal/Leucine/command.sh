python=$1
shift
params=$*
$python ../../../pdielec $params -matrix ptfe -sigma 5 \
                 -plate 0 0 1 -ellipsoid 0 0 1 0.5 -plate 1 0 0 \
                 -optical_tensor 2.3471 0 0.0786 0 2.3942 0 0.0786 0 2.2621 \
                 -csv command.csv \
                 -masses program LEUCINE_FREQUENCY_PBED3_631Gdp_FULLOPTIMIZATON.out  $*
