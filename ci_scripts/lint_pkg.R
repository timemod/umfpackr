library(lintr)

lints <- lint_package("pkg")

if (length(lints) > 0) {
  cat("\nLinter detected style issues:\n")
  print(lints)
  stop("Linter failed. Please fix the above issues.")
} else {
  message("No linting issues found")
}
