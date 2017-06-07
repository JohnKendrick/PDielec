python ../../../pdielec -matrix ptfe -sigma 5 -method maxwell -method bruggeman -method ap -vf 0.1\
                         -needle 0 0 1 -ellipsoid 0 0 1 0.5 -plate 1 0 0 \
                         -LO 1 1 1 \
                         -csv command1.csv \
                         -neutral -program vasp OUTCAR  $*
python ../../../pdielec -matrix ptfe -sigma 5 -method maxwell -vf 0.1\
                         -csv command2.csv -masses average -print \
                         -eckart -program vasp OUTCAR  $*
cat command1.csv command2.csv > command.csv
rm -f command1.csv command2.csv
