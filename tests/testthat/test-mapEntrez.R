test_that("mapEntrez works", {
  expect_identical(mapEntrez(c("A1BG", "A1CF")), c("A1BG", "A1CF"))
  expect_identical(mapEntrez(c("1", "29974")), c("A1BG", "A1CF"))
})

test_that("mapHugo works", {
  expect_identical(mapHugo(c("A1BG", "A1CF")), c("1", "29974"))
  expect_identical(mapHugo(c("1", "29974")), c("1", "29974"))
})
