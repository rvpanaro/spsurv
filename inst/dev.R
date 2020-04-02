setwd("~/Documents/spsurv")

devtools::load_all()
devtools::document()
# devtools::install()
devtools::install(quick=TRUE)

usethis::use_build_ignore("_pkgdown.yml")
usethis::use_build_ignore("codecov.yml")
usethis::use_build_ignore(".travis.yml")
usethis::use_build_ignore("_config.yml")
usethis::use_build_ignore("cran-comments.md")
usethis::use_build_ignore("inst/dev.R")
usethis::use_build_ignore("inst/load.R")

#devtools::build_vignettes()
devtools::check()
devtools::check_man()
devtools::missing_s3()
devtools::release_checks()
# devtools::check_win_devel()

# usethis::use_revdep()
# devtools::revdep()

devtools::spell_check()
devtools::test()

# devtools::reload()
t1 <- Sys.time(); devtools::run_examples();
t2 <- Sys.time(); t2-t1

devtools::build_manual()
covr::codecov(token = "5ac8c633-0916-4cdc-9b1f-7391b8dafe7f")

library("hexSticker")
library(showtext)
font_add_google("Russo One", "f")
sticker("inst/figures/Rplot.png", package="spsurv", p_size=20, s_x=1, s_y=.75, s_width=.6,
        h_fill="#01ACB6", h_color="#9be203",
        filename="inst/figures/imgfile.png", p_family = "f")
devtools::submit_cran()
