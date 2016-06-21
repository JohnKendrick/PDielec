python ../../../pdielec -matrix ptfe -vmin 900 -vmax 1100 -sigma 2 -sphere -plate 1 0 0 -plate -1 1 0 -needle 0 0 1 -ellipsoid 0 0 1 2 -method ap -method maxwell -csv command.csv .  $*
