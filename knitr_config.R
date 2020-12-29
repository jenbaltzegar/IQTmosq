#options(knitr.kable.NA='')
knitr::opts_chunk$set(
    ## over-ride latex page layout / float positioning
    fig.pos='h!',
    warning=FALSE,
    message=FALSE,
    eval=TRUE,
    include=TRUE,
    echo=FALSE,
    fig.width=10,
    #out.width='12in',
    #fig.height=4,
    dpi=300,
    dev=c('pdf','png'),
    ## beamer text (not captions)
    #size = 'tiny',
    size = 'footnotesize'
)
