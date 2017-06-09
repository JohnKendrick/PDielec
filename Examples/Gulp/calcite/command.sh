python ../../../pdielec -matrix ptfe -method ap -method maxwell -sphere -plate -1 -1 -2 -vmin 0 -vmax 2000 -vf 0.1 calcite.gout  -masses program -csv command1.csv  $*
python ../../../pdielec -matrix ptfe -method ap -method maxwell -sphere -vf 0.1 calcite.gout -masses isotope -print  -csv command2.csv  $*
cat command1.csv command2.csv > command.csv
rm -f command1.csv command2.csv
