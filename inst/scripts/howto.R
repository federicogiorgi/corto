### Git tricks
# git config --global user.name "federicogiorgi"
# git config --global user.email "federico.giorgi@gmail.com"

### Switching remote URLs from SSH to HTTPS
# # git remote -v
# >origin  https://github.com/federicogiorgi/corto (fetch)
# >origin  https://github.com/federicogiorgi/corto (push)
# git remote set-url origin https://github.com/federicogiorgi/corto.git
### Switching remote URLS from HTTPS to SSH
# git remote set-url origin git@github.com:federicogiorgi/corto.git


# if (!require(git2r)) {
#   install.packages("git2r", repos="http://cran.r-project.org/", quiet=TRUE)
# }
# library(git2r)
# repo<-repository(".")
# push(repo, credentials=cred_user_pass("federicogiorgi",pass))


# Building manuals and namespace using roxygen2
library(devtools)
unlink("NAMESPACE")
document()
cat('importFrom("stats", "cor", "pt", "qt", "setNames", "var",
    "pchisq","ks.test", "pnorm", "qnorm", "quantile")\n',file="NAMESPACE",append=TRUE)
cat('importFrom("utils", "setTxtProgressBar", "txtProgressBar")\n',file="NAMESPACE",append=TRUE)
cat('import("doSNOW")\n',file="NAMESPACE",append=TRUE)
cat('import("snow")\n',file="NAMESPACE",append=TRUE)
cat('import("foreach")\n',file="NAMESPACE",append=TRUE)
cat('importFrom("dplyr", "%>%", "filter", "group_by")\n',file="NAMESPACE",append=TRUE)
cat('importFrom("graphics", "abline", "grid", "layout", "legend", "lines",
           "par", "plot", "rect", "text")\n',file="NAMESPACE",append=TRUE)
build()
check(document=FALSE)

# # Remove .Rhistory
# unlink(".Rhistory")
# unlink("R/.Rhistory")

