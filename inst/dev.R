setwd("~/Documents/spsurv")

devtools::load_all()
devtools::document()
# devtools::install()
devtools::install(quick=TRUE)

usethis::use_build_ignore(".travis.yml")
usethis::use_build_ignore("_config.yml")
usethis::use_build_ignore("cran-comments.md")
usethis::use_build_ignore("inst/dev.R")

#devtools::build_vignettes()
devtools::check()
devtools::check_man()
devtools::missing_s3()
devtools::release_checks()
devtools::check_win_devel()

usethis::use_revdep()
devtools::revdep()

devtools::spell_check()
devtools::test()

# devtools::reload()
devtools::run_examples()
devtools::build_manual()
devtools::submit_cran()
#install.packages("crossSurv_0.0.0.9000.tar.gz", repos = FALSE, dependencies=TRUE)
# para eliminar o sequinte erro: Found the following hidden files and directories: .travis.yml
# git add .travis.yml .Rbuildignore
# git commit -m 'enable continuous integration via craigcitro/r-travis'
# git push


# for citation:
# citation(auto = packageDescription("YPPE"))
