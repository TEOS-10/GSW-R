# install.packages("codemetar")
requireNamespace(c("codemetar", "devtools", "urlchecker", "rhub", "revdepcheck"))
# codemeta changes a timestamp, so requiring a commit after every call. That is
# senseless, so I only run the false part of the following conditional in the
# run-up to a release.
if (FALSE) {
    codemetar::write_codemeta()
} else {
    message("run 'codemetar::write_codemeta()' and then git push")
}
t <- devtools::spell_check()
stopifnot(t == "No spelling errors found.")
urlchecker::url_check()
devtools::check_win_release()
devtools::check_win_devel()
devtools::check_win_oldrelease()

# 2024-08-19 DEK: rhub has been replaced by rhub2, so the following no longer
# work. The new system is complicated, and has the disadvantage that it only
# works on code that has been pushed to github ... so now we rely on github
# actions.
#    rhub::check_for_cran(email = "Dan.Kelley@Dal.Ca", show_status = FALSE)
#    rhub::check(platform = "debian-clang-devel", show_status = FALSE)
remotes::install_github("r-lib/revdepcheck")
revdepcheck::revdep_reset()
revdepcheck::revdep_check(num_workers = 4)
