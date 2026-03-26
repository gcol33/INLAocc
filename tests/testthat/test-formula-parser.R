# =============================================================================
# Formula parser tests — lme4-style random effects syntax
#
# Tests the AST-walking parser (findbars_ast, nobars_ast, bar_to_occu_re,
# parse_re_formula) for all supported lme4 random effects patterns.
# No INLA dependency needed — these are pure formula-parsing tests.
# =============================================================================

# --- Helper: shortcut to parse and return re_list ---
parse_re <- function(formula) {
  parse_re_formula(formula)$re_list
}

parse_fixed <- function(formula) {

  parse_re_formula(formula)$fixed
}

# ==========================================================================
# Basic features (regression)
# ==========================================================================

test_that("(1 | group) produces a single intercept RE", {
  re <- parse_re(~ x + (1 | group))
  expect_length(re, 1)
  expect_s3_class(re[[1]], "occu_re")
  expect_equal(re[[1]]$type, "intercept")
  expect_equal(re[[1]]$group, "group")
  expect_true(re[[1]]$correlated)
})

test_that("(1 + x | group) produces intercept + slope RE", {
  re <- parse_re(~ (1 + x | group))
  expect_length(re, 2)
  expect_equal(re[[1]]$type, "intercept")
  expect_equal(re[[1]]$group, "group")
  expect_equal(re[[2]]$type, "slope")
  expect_equal(re[[2]]$group, "group")
  expect_equal(re[[2]]$covariate, "x")
  expect_true(re[[1]]$correlated)
  expect_true(re[[2]]$correlated)
})

test_that("(1 | a/b) produces nested grouping a:b", {
  re <- parse_re(~ (1 | a/b))
  expect_length(re, 1)
  expect_equal(re[[1]]$group, "a:b")
  expect_equal(re[[1]]$type, "intercept")
})

test_that("crossed random effects (1|g1) + (1|g2) produce two RE objects", {
  re <- parse_re(~ x + (1 | g1) + (1 | g2))
  expect_length(re, 2)
  expect_equal(re[[1]]$group, "g1")
  expect_equal(re[[2]]$group, "g2")
  expect_equal(re[[1]]$type, "intercept")
  expect_equal(re[[2]]$type, "intercept")
})

test_that("fixed-effects formula is correctly stripped of bar terms", {
  fixed <- parse_fixed(~ x + z + (1 | group))
  expect_equal(sort(attr(terms(fixed), "term.labels")), c("x", "z"))
})

# ==========================================================================
# Intercept suppression
# ==========================================================================

test_that("(0 + x | group) suppresses intercept, keeps slope", {
  re <- parse_re(~ (0 + x | group))
  expect_length(re, 1)
  expect_equal(re[[1]]$type, "slope")
  expect_equal(re[[1]]$covariate, "x")
  expect_equal(re[[1]]$group, "group")
})

test_that("(-1 + x | group) suppresses intercept via -1", {
  re <- parse_re(~ (-1 + x | group))
  expect_length(re, 1)
  expect_equal(re[[1]]$type, "slope")
  expect_equal(re[[1]]$covariate, "x")
  expect_equal(re[[1]]$group, "group")
})

test_that("(0 + x1 + x2 | group) gives two slopes, no intercept", {
  re <- parse_re(~ (0 + x1 + x2 | group))
  expect_length(re, 2)
  expect_equal(re[[1]]$type, "slope")
  expect_equal(re[[1]]$covariate, "x1")
  expect_equal(re[[2]]$type, "slope")
  expect_equal(re[[2]]$covariate, "x2")
})

# ==========================================================================
# Interaction grouping
# ==========================================================================

test_that("(1 | group1:group2) uses interaction as grouping variable", {
  re <- parse_re(~ (1 | group1:group2))
  expect_length(re, 1)
  expect_equal(re[[1]]$type, "intercept")
  # The parser deparses the RHS; interaction expressions become "group1:group2"
  expect_equal(re[[1]]$group, "group1:group2")
})

# ==========================================================================
# Uncorrelated random effects (|| double-bar syntax)
# ==========================================================================

test_that("(1 || group) sets correlated = FALSE", {
  re <- parse_re(~ (1 || group))
  expect_length(re, 1)
  expect_equal(re[[1]]$type, "intercept")
  expect_equal(re[[1]]$group, "group")
  expect_false(re[[1]]$correlated)
})

test_that("(1 + x || group) sets correlated = FALSE on both RE objects", {
  re <- parse_re(~ (1 + x || group))
  expect_length(re, 2)
  expect_equal(re[[1]]$type, "intercept")
  expect_equal(re[[2]]$type, "slope")
  expect_equal(re[[2]]$covariate, "x")
  expect_false(re[[1]]$correlated)
  expect_false(re[[2]]$correlated)
})

test_that("|| is stripped from fixed formula like |", {
  parsed <- parse_re_formula(~ x + (1 || group))
  fixed <- parsed$fixed
  terms_labels <- attr(terms(fixed), "term.labels")
  expect_equal(terms_labels, "x")
  expect_length(parsed$re_list, 1)
})

# ==========================================================================
# Implicit intercept (bare covariate on LHS)
# ==========================================================================

test_that("(x | group) treats x as a slope term with implicit intercept", {
  re <- parse_re(~ (x | group))
  # lme4 convention: bare variable on LHS is a slope; intercept is implicit
  expect_length(re, 2)
  expect_equal(re[[1]]$type, "intercept")
  expect_equal(re[[2]]$type, "slope")
  expect_equal(re[[2]]$covariate, "x")
  expect_equal(re[[1]]$group, "group")
  expect_equal(re[[2]]$group, "group")
})

# ==========================================================================
# occu_re() constructor — correlated parameter
# ==========================================================================

test_that("occu_re() defaults correlated to TRUE", {
  re <- occu_re("intercept", group = "g")
  expect_true(re$correlated)
})

test_that("occu_re() accepts correlated = FALSE", {
  re <- occu_re("intercept", group = "g", correlated = FALSE)
  expect_false(re$correlated)
})

# ==========================================================================
# Edge cases
# ==========================================================================

test_that("formula with no random effects returns empty re_list", {
  parsed <- parse_re_formula(~ x + z)
  expect_length(parsed$re_list, 0)
  expect_equal(sort(attr(terms(parsed$fixed), "term.labels")), c("x", "z"))
})

test_that("intercept-only formula with RE works", {
  parsed <- parse_re_formula(~ (1 | group))
  expect_length(parsed$re_list, 1)
  # Fixed formula should be ~ 1
  expect_length(attr(terms(parsed$fixed), "term.labels"), 0)
})

test_that("mixed | and || in same formula are parsed independently", {
  re <- parse_re(~ (1 | g1) + (1 + x || g2))
  expect_length(re, 3)
  # First: correlated intercept on g1

  expect_equal(re[[1]]$group, "g1")
  expect_true(re[[1]]$correlated)
  # Second: uncorrelated intercept on g2
  expect_equal(re[[2]]$group, "g2")
  expect_false(re[[2]]$correlated)
  # Third: uncorrelated slope on g2
  expect_equal(re[[3]]$group, "g2")
  expect_equal(re[[3]]$covariate, "x")
  expect_false(re[[3]]$correlated)
})
