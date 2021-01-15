library(knitr)
## depends on / assumes:
# source('make_jb.R')
knit('doc.Rnw')
tinytex::pdflatex('doc.tex', pdf_file='output/doc.pdf')
