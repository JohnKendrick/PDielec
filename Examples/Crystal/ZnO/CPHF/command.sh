python ../../../../pdielec -vmin 300 -vmax 800 -matrix ptfe -sigma 5 \
                            -sphere -needle 0 0 1 -plate 1 0 0 \
                            -csv command.csv \
                            -hessian crystal -program crystal ZnO_CPHF.out   $*
