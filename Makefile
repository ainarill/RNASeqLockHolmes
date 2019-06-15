all : index.html

index.html : index.Rmd QAanalysis.knit.md DEanalysis.knit.md FEanalysis.knit.md QAanalysisTumor.knit.md DEanalysisTumor.knit.md FEanalysisTumor.knit.md
	Rscript -e "rmarkdown::render('$<')"

QAanalysis.knit.md : QAanalysis.Rmd
	Rscript -e "rmarkdown::render('$<', run_pandoc=FALSE, clean=FALSE)"

DEanalysis.knit.md : DEanalysis.Rmd
	Rscript -e "rmarkdown::render('$<', run_pandoc=FALSE, clean=FALSE)"

FEanalysis.knit.md : FEanalysis.Rmd
	Rscript -e "rmarkdown::render('$<', run_pandoc=FALSE, clean=FALSE)"

QAanalysisTumor.knit.md : QAanalysisTumor.Rmd
	Rscript -e "rmarkdown::render('$<', run_pandoc=FALSE, clean=FALSE)"

DEanalysisTumor.knit.md : DEanalysisTumor.Rmd
	Rscript -e "rmarkdown::render('$<', run_pandoc=FALSE, clean=FALSE)"

FEanalysisTumor.knit.md : FEanalysisTumor.Rmd
	Rscript -e "rmarkdown::render('$<', run_pandoc=FALSE, clean=FALSE)"