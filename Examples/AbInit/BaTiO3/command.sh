python ../../../pdielec -method ap -method bruggeman -matrix ptfe -sigma 5 \
                 -plate 0 0 1 -ellipsoid 0 0 1 0.5 -plate 1 0 0 \
                 -LO 1 0 0 -LO 0 1 0 -LO 0 0 1 \
                 -csv command.csv \
                 -program abinit BaTiO3.out $*
