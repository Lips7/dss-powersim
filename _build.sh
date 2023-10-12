Rscript -e "bookdown::render_book(output_format='bookdown::gitbook')"
cp .nojekyll docs/
rm -rf *.rds *.log