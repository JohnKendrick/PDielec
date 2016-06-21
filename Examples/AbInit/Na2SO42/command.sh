python ../../../pdielec -method ap -method maxwell -matrix ptfe -sigma 5 \
                 -plate 0 0 1 -ellipsoid 0 0 1 0.5 -plate 1 0 0 \
                 -LO_cart 1 0 0 -LO_cart 0 1 0 -LO_cart 0 0 1 \
                 -LO 1 0 0 -LO 0 1 0 -LO 0 0 1 \
                 -csv command.csv \
                 -program abinit Na2SO42.out $*
