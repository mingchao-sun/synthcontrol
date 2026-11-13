# Panel validation and messages for bad input
# Converted from tests/test_stability.R. Runs from the module folder regardless of the working directory.
root <- normalizePath(file.path(testthat::test_path(), "..", ".."))

testthat::test_that("Panel validation and messages for bad input", {
  withr::local_dir(root)
  ok <- function(cond, msg) testthat::expect_true(isTRUE(cond), info = msg)
  check <- function(name, cond) testthat::expect_true(isTRUE(cond), info = name)
  # Stability tests for the panel validation (scm_build_panel).
  # These mirror what the jamovi analysis does with unusual input: missing values,
  # text in the time variable, unusual unit and variable names, empty data,
  # duplicate rows, and invalid treatment settings. Each bad case must stop with
  # a clear message and must not crash.
  source("R/engine.R")

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

  cat("16. Missing unit and time values get their own messages\n")
  du <- d; du$state[3] <- NA
  check("blank unit is named", grepl("unit variable has missing", err_msg(run(du))))
  dt <- d; dt$year[5] <- NA
  check("blank time is named", grepl("time variable has missing", err_msg(run(dt))))
  dtx <- d; dtx$year <- as.character(dtx$year); dtx$year[2] <- "19x0"
  check("non-numeric time is named", grepl("must be numeric", err_msg(run(dtx))))

  cat("17. Predictor with no pre-treatment values\n")
  dp <- d; dp$lnincome <- 1
  dp$lnincome[dp$state == "Utah" & dp$year <= 1988] <- NA
  Xp_msg <- err_msg(scm_predictor_matrix(dp, "state", "year", "lnincome", unique(dp$state),
                                         1970:1988, stat = "mean"))
  check("all-missing predictor names the predictor", grepl("Predictor 'lnincome'", Xp_msg))
  check("message names the unit", grepl("Utah", Xp_msg))
  dl <- d; dl$lnincome <- 1
  dl$lnincome[dl$state == "Utah" & dl$year == 1988] <- NA
  Xl_msg <- err_msg(scm_predictor_matrix(dl, "state", "year", "lnincome", unique(dl$state),
                                         1970:1988, stat = "last"))
  check("missing last pre-treatment value is caught", grepl("Utah", Xl_msg))
  Xok <- scm_predictor_matrix(dl, "state", "year", "lnincome", unique(dl$state), 1970:1987, stat = "last")
  check("last value from an earlier window is fine", is.matrix(Xok) && all(is.finite(Xok)))

  cat("18. Scenario offsets are checked against the data\n")
  times <- 1970:2000
  check("valid scenarios pass", isTRUE(scm_check_scenarios(times, 1989, c(-2, 1), c(0, 2), c(-6, -2))))
  m <- err_msg(scm_check_scenarios(times, 1989, intime = 6))
  check("positive in-time offset is rejected with its value", grepl("offset 6", m))
  m <- err_msg(scm_check_scenarios(times, 1989, timing = -25))
  check("timing offset leaving too few periods is rejected", grepl("Treatment timing offset -25", m))
  m <- err_msg(scm_check_scenarios(times, 1989, timing = 20))
  check("timing offset past the data is rejected", grepl("Treatment timing offset 20", m))
  m <- err_msg(scm_check_scenarios(times, 1989, drops = 18))
  check("pre-period drop leaving fewer than two periods is rejected", grepl("Pre-period drop 18", m))
  m <- err_msg(scm_check_scenarios(times, 1989, drops = 1.5))
  check("non-integer drop is rejected", grepl("whole numbers", m))
  check("parse drops duplicates", identical(scm_parse_list("2, 2, 4", "x"), c(2, 4)))
  check("fractional labels format without error", grepl("1988.5", scm_lab_intime(1988.5)))

  cat("19. Placebo checkpoint is called once per donor\n")
  Yp <- matrix(0, 1, 1)
  Ys <- as.matrix(d[d$year <= 1999, ] |> transform(cs = ave(cigsale, state, FUN = identity)) |> with(0))
  Ym <- { w <- reshape(d[, c("year", "state", "cigsale")], idvar = "year", timevar = "state", direction = "wide")
          m <- as.matrix(w[, -1]); rownames(m) <- w$year; colnames(m) <- sub("cigsale.", "", colnames(m)); m }
  pre_i <- which(as.integer(rownames(Ym)) <= 1988); post_i <- which(as.integer(rownames(Ym)) >= 1989)
  calls <- 0
  pl_ck <- scm_placebo(Ym, "California", pre_i, post_i, starts = 1,
                       checkpoint = function() calls <<- calls + 1)
  check("checkpoint called once per placebo donor", calls == 38)
  check("placebo results have one entry per donor", length(pl_ck$placebo) == 38)
  slim <- scm_placebo_slim(pl_ck)
  check("slim placebo gives the same p-value", identical(scm_pvalue(slim), scm_pvalue(pl_ck)))

  cat("20. Loo subset and timing rows match the table labels\n")
  tm <- scm_timing(Ym, "California", 1989, offsets = c(-2, 1), starts = 1)
  check("timing labels match the .init row labels",
        identical(tm$scenario, c(scm_lab_timing(1987), scm_lab_timing(1990))))

})
