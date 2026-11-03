# BayesPrism on the prepared inputs: cell.type = major, cell.state = fine; theta from the initial Gibbs pass (type and state)
suppressMessages(library(BayesPrism))
args <- commandArgs(trailingOnly = TRUE); name <- args[1]; O <- "/home/claude/work/prop/bp"
ref <- as.matrix(read.csv(gzfile(file.path(O, paste0(name, "_ref.csv.gz"))), row.names = 1, check.names = FALSE))
lab <- read.csv(file.path(O, paste0(name, "_ref_labels.csv")))
mix <- as.matrix(read.csv(gzfile(file.path(O, paste0(name, "_mix.csv.gz"))), row.names = 1, check.names = FALSE))
storage.mode(ref) <- "numeric"; storage.mode(mix) <- "numeric"
pr <- new.prism(reference = ref, mixture = mix, input.type = "count.matrix", cell.type.labels = lab$major,
                cell.state.labels = lab$fine, key = NULL, outlier.cut = 0.01, outlier.fraction = 0.1)
bp <- run.prism(prism = pr, n.cores = 2, update.gibbs = FALSE)   # no tumour key: one Gibbs pass, theta "first"
write.csv(get.fraction(bp, which.theta = "first", state.or.type = "type"), file.path(O, paste0(name, "_theta_type.csv")))
write.csv(get.fraction(bp, which.theta = "first", state.or.type = "state"), file.path(O, paste0(name, "_theta_state.csv")))
cat("done", name, "\n")
