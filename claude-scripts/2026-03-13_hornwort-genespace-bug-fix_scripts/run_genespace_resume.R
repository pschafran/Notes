library(GENESPACE)
library(data.table)

wd <- "/media/data/projects/hornwort_sex_chromosomes/analysis/synteny/genespace_hornworts_20260312/"

# Apply patch for zero-block MCScanX crash (GENESPACE v1.3.1 bug)
source(file.path(wd, "patch_run_mcscanx.R"))

# Resume from saved session
gpar <- init_genespace(
    wd = wd,
    path2mcscanx = "~/bin/MCScanX/",
    nCores = 1
)

results <- run_genespace(gsParam = gpar)
