args <- commandArgs(trailingOnly = TRUE)
output_html <- if (length(args) >= 1L) args[[1]] else "coverage.html"
output_xml <- if (length(args) >= 2L) args[[2]] else "coverage.xml"

escape_html <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}

html_table <- function(data) {
  headers <- paste(sprintf("<th>%s</th>", escape_html(names(data))), collapse = "")
  rows <- apply(data, 1, function(row) {
    cells <- paste(sprintf("<td>%s</td>", escape_html(as.character(row))), collapse = "")
    sprintf("<tr>%s</tr>", cells)
  })
  paste0(
    "<table><thead><tr>", headers, "</tr></thead><tbody>",
    paste(rows, collapse = ""),
    "</tbody></table>"
  )
}

cov <- covr::package_coverage()
print(cov)
cat(as.character(covr::to_cobertura(cov)), file = output_xml)

coverage_data <- as.data.frame(cov)
coverage_data$covered <- coverage_data$value > 0

file_groups <- split(coverage_data, coverage_data$filename)
file_summary <- do.call(rbind, lapply(file_groups, function(group) {
  data.frame(
    file = group$filename[[1]],
    lines = nrow(group),
    covered = sum(group$covered),
    uncovered = sum(!group$covered),
    coverage = sprintf("%.2f%%", 100 * mean(group$covered)),
    stringsAsFactors = FALSE
  )
}))
row.names(file_summary) <- NULL
file_summary <- file_summary[order(file_summary$uncovered, file_summary$file, decreasing = TRUE), ]

uncovered_lines <- unique(coverage_data[!coverage_data$covered, c("filename", "functions", "first_line")])
names(uncovered_lines) <- c("file", "function", "line")
uncovered_lines[["function"]][uncovered_lines[["function"]] == ""] <- "<top-level>"
uncovered_lines <- uncovered_lines[order(uncovered_lines$file, uncovered_lines$line), ]

overall_coverage <- sprintf("%.2f%%", covr::percent_coverage(cov))
generated_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")

html <- paste0(
  "<!DOCTYPE html><html><head><meta charset='utf-8'>",
  "<title>kdensity Coverage</title>",
  "<style>",
  "body{font-family:system-ui,-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;margin:2rem;color:#17212b;background:#f6f8fb;}",
  "main{max-width:1100px;margin:0 auto;background:#fff;padding:2rem 2.5rem;border-radius:16px;box-shadow:0 12px 30px rgba(23,33,43,.08);}",
  "h1,h2{margin:0 0 1rem;}p{line-height:1.5;}table{border-collapse:collapse;width:100%;margin:1rem 0 2rem;}",
  "th,td{padding:.65rem .8rem;border-bottom:1px solid #d8e0ea;text-align:left;vertical-align:top;}",
  "th{background:#eef3f8;font-weight:600;}.meta{display:flex;gap:2rem;flex-wrap:wrap;margin:1rem 0 2rem;}",
  ".card{background:#eef3f8;padding:1rem 1.25rem;border-radius:12px;}",
  "</style></head><body><main>",
  "<h1>kdensity Coverage Report</h1>",
  sprintf("<p>Generated %s.</p>", escape_html(generated_at)),
  "<div class='meta'>",
  sprintf("<div class='card'><strong>Overall coverage</strong><br>%s</div>", escape_html(overall_coverage)),
  sprintf("<div class='card'><strong>Files</strong><br>%s</div>", nrow(file_summary)),
  sprintf("<div class='card'><strong>Uncovered lines</strong><br>%s</div>", nrow(uncovered_lines)),
  "</div>",
  "<h2>Coverage by File</h2>",
  html_table(file_summary),
  "<h2>Uncovered Lines</h2>",
  if (nrow(uncovered_lines) > 0) html_table(uncovered_lines) else "<p>All tracked lines are covered.</p>",
  "</main></body></html>"
)

writeLines(html, con = output_html)
