test_that("Variable regions are detected", {
  expect_length(summarize_variable_regions(unitorders_example, 0.75), 2)
})

test_that("Consensus sequence output is expected length",{
  expect_length(call_consensus_sequence(unitorders_example, unit_frequency_example, unit2ascii),
                2)
})
