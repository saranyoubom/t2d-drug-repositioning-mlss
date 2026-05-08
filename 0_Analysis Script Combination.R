# compile_scripts.R
# Compiles all R scripts in the working directory into a single output file

compile_r_scripts <- function(
    output_file = "compiled_output.R",
    working_dir = getwd(),
    exclude_self = TRUE,
    add_separators = TRUE
) {
  # Get all .R files in the working directory
  r_files <- list.files(
    path = working_dir,
    pattern = "\\.R$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  # Exclude this script itself (optional)
  if (exclude_self) {
    self_name <- basename(output_file)
    r_files <- r_files[basename(r_files) != "compile_scripts.R"]
    r_files <- r_files[basename(r_files) != self_name]
  }
  
  if (length(r_files) == 0) {
    message("No R files found in: ", working_dir)
    return(invisible(NULL))
  }
  
  message("Found ", length(r_files), " R file(s) to compile:")
  for (f in r_files) message("  - ", basename(f))
  
  # Open output connection
  output_path <- file.path(working_dir, output_file)
  out_con <- file(output_path, open = "w")
  on.exit(close(out_con))
  
  # Write header
  writeLines(c(
    "# ============================================================",
    paste0("# Compiled R Script"),
    paste0("# Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("# Source directory: ", working_dir),
    paste0("# Files compiled: ", length(r_files)),
    "# ============================================================",
    ""
  ), con = out_con)
  
  # Append each file's content
  for (f in r_files) {
    if (add_separators) {
      writeLines(c(
        "",
        "# ============================================================",
        paste0("# Source file: ", basename(f)),
        "# ============================================================",
        ""
      ), con = out_con)
    }
    
    # Read and write the file contents
    lines <- readLines(f, warn = FALSE)
    writeLines(lines, con = out_con)
    writeLines("", con = out_con)  # trailing newline between files
  }
  
  message("\nCompiled successfully -> ", output_path)
  invisible(output_path)
}

# Run it
compile_r_scripts(
  output_file = "compiled_output.R",
  working_dir = getwd(),
  exclude_self = TRUE,
  add_separators = TRUE
)

