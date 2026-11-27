for (f in c("mlss_port.R", "wp3_mlss.R", "wp1_ins1_analysis.R")) {
  r <- tryCatch({ parse(file = file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), f)); "ok" }, error = function(e) conditionMessage(e))
  cat(sprintf("%-22s %s\n", f, r))
}
