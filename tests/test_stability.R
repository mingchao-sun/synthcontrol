# Stability tests for the panel validation (scm_build_panel).
# These mirror what the jamovi analysis does with unusual input: missing values,
# text in the time variable, unusual unit and variable names, empty data,
# duplicate rows, and invalid treatment settings. Each bad case must stop with
# a clear message and must not crash.
source("R/engine.R")

ok <- 0; fail <- 0
check <- function(name, cond) {
    if (isTRUE(cond)) { ok <<- ok + 1; cat("  ok -", name, "\n") }
    else { fail <<- fail + 1; cat("  FAIL -", name, "\n") }
}
err_msg <- function(expr) tryCatch({ expr; NA_character_ }, error = function(e) conditionMessage(e))

d <- read.csv("data/prop99_long.csv", stringsAsFactors = FALSE)
base_args <- list(Y_name = "cigsale", unit_name = "state", time_name = "year",
                  treated = "California", treatmentTime = 1989)
run <- function(df, treated = base_args$treated, exclude = character(0), ...) {
  a <- modifyList(base_args, list(...))
  a$treated <- treated; a$exclude <- exclude
  do.call(scm_build_panel, c(list(df = df), a))
}

cat("1. Valid panel\n")
p <- run(d)
check("Prop 99 panel builds", is.list(p))
check("19 pre and 12 post periods", length(p$pre) == 19 && length(p$post) == 12)
check("39 units", ncol(p$Y) == 39)

cat("2. Missing outcome\n")
d_na <- d; d_na$cigsale[5] <- NA
m <- err_msg(run(d_na))
check("missing outcome is rejected with a balance message", grepl("balanced", m))

cat("3. Non-numeric time\n")
d_txt <- d; d_txt$year <- paste0("y", d_txt$year)
m <- err_msg(run(d_txt))
check("non-numeric time is rejected with a numeric message", grepl("numeric", m))

cat("4. Time stored as text but numeric values\n")
d_chr <- d; d_chr$year <- as.character(d_chr$year)
check("character years are coerced and accepted", is.list(run(d_chr)))

cat("5. Unusual unit names\n")
d_sp <- d; d_sp$state <- gsub("California", "New York", d_sp$state)
d_sp$state <- gsub("Utah", "A&B (Utah)", d_sp$state)
p2 <- run(d_sp, treated = "New York")
check("unit names with spaces and symbols work", "A&B (Utah)" %in% colnames(p2$Y))

cat("6. Unusual outcome column name\n")
d_nm <- d; names(d_nm)[names(d_nm) == "cigsale"] <- "cig sales"
p3 <- run(d_nm, Y_name = "cig sales")
check("outcome name with a space works", is.list(p3))

cat("7. Empty data\n")
m <- err_msg(run(d[0, ]))
check("empty data is rejected without crashing", !is.na(m))

cat("8. Duplicate unit-time rows\n")
d_dup <- rbind(d, d[1, ])
m <- err_msg(run(d_dup))
check("duplicate unit-time rows are rejected", grepl("duplicate", m))

cat("9. Treated unit not in data\n")
m <- err_msg(run(d, treated = "Texas_missing"))
check("unknown treated unit is rejected", grepl("treated unit", m))

cat("10. Treated unit in donor exclusion\n")
m <- err_msg(run(d, exclude = "California"))
check("treated unit cannot be excluded", grepl("cannot appear", m))

cat("11. Unknown unit in donor exclusion\n")
m <- err_msg(run(d, exclude = "Atlantis"))
check("unknown exclusion is rejected", grepl("Unknown unit", m))

cat("12. Treatment too early\n")
m <- err_msg(run(d, treatmentTime = 1971))
check("fewer than two pre-periods is rejected", grepl("pre-treatment", m))

cat("13. Treatment after the last period\n")
m <- err_msg(run(d, treatmentTime = 2050))
check("no post-period is rejected", grepl("post-treatment", m))

cat("14. Too few donors after exclusion\n")
m <- err_msg(run(d, exclude = setdiff(unique(d$state), c("California", "Utah"))))
check("fewer than two donors is rejected", grepl("donor", m))

cat("15. Donor exclusion works\n")
p4 <- run(d, exclude = c("Utah", "Nevada"))
check("excluded donors are dropped", !("Utah" %in% p4$units) && ncol(p4$Y) == 37)

cat(sprintf("\n%d passed, %d failed\n", ok, fail))
if (fail > 0) quit(status = 1)
cat("All stability tests passed.\n")
