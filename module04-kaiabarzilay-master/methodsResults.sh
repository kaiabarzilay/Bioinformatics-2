#!/usr/bin/env bash
R -e "rmarkdown::render('methodsResults.html', output_format='all')"
