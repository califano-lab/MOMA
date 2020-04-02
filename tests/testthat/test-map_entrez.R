test_that("map_entrez works", {
  expect_identical(map_entrez(c("A1BG", "A1CF")), c("A1BG", "A1CF"))
  expect_identical(map_entrez(c("1", "29974")), c("A1BG", "A1CF"))
})

test_that("map_hugo works", {
  expect_identical(map_hugo(c("A1BG", "A1CF")), c("1", "29974"))
  expect_identical(map_hugo(c("1", "29974")), c("1", "29974"))
})
