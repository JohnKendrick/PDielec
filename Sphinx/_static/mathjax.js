MathJax.Hub.Config({
    TeX: {
        extensions: ['AMSmath.js'] ,
        Macros: { 
            bm:       ['{\\boldsymbol{#1}}', 1],
            tensor:   ['{\\bar{\\bar{#1}}}', 1],
            tensorbs: ['{\\bar{\\bar{\\bm{#1}}}}', 1],
            tensorbf: ['{\\bar{\\bar{\\bm{#1}}}}', 1],
            fieldbf:  ['{\\bar{\\bm{#1}}}', 1],
            water: "\\H_{2}O"
        }
    }
});
