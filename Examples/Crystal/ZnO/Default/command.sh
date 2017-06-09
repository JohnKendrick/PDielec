python ../../../../pdielec -vmin 300 -vmax 800 -matrix ptfe -sigma 5 \
                            -sphere -needle 0 0 1 -plate 1 0 0 \
                            -optical 5.73897 5.73897 5.67346  -csv command.csv \
                            -masses program -hessian crystal -program crystal ZnO_default.out   $*
