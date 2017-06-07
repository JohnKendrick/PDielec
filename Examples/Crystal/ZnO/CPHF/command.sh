python ../../../../pdielec -vmin 300 -vmax 800 -matrix ptfe -sigma 5 \
                            -sphere -needle 0 0 1 -plate 1 0 0 \
                            -csv command1.csv \
                            -hessian crystal -program crystal ZnO_CPHF.out   $*
python ../../../../pdielec -vmin 300 -vmax 800 -matrix ptfe -sigma 5 \
                            -sphere \
                            -csv command2.csv -masses program -print \
                            -program crystal ZnO_CPHF.out   $*
cat command1.csv command2.csv > command.csv
rm -f command1.csv command2.csv
