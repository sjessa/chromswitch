context("External cluster validity functions")


test_that("External cluster validity indices are correct", {

    metadata <- data.frame(Sample = c("A", "B", "C", "D"),
                           Group = c(1, 1, 2, 2),
                           stringsAsFactors = FALSE)

    clusters <- c(A = 1, B = 1, C = 0, D = 0)

    expect_equal(as.data.frame(contingencyTable(clusters, metadata)),
                data.frame(Var1 = c("1", "2", "1", "2"),
                Var2 = c("0", "0", "1", "1"),
                Freq = c(0, 2, 2, 0)))


    clusters <- c(0, 0, 2, 1, 1, 0, 1)
    classes <- c("A", "A", "A", "B", "B", "A", "B")
    ct <- table(classes, clusters)

    expect_equal(purity(contingency = ct), 1)
    expect_equal(purity(c = classes, k = clusters), 1)

    expect_equal(classEntropy(ct), 0.6829081)
    expect_equal(clusterEntropy(ct), 1.004242, tolerance = 1e-6)

    expect_equal(conditionalClassEntropy(ct), 0)
    expect_equal(conditionalClusterEntropy(ct), 0.3213344, tolerance = 1e-6)

    expect_equal(homogeneity(ct), 1)
    expect_equal(completeness(ct), 0.6800231, tolerance = 1e-6)
    expect_equal(vMeasure(ct), 0.8095402)

    expect_equal(NMI(clusters, classes), 0.8246351, tolerance = 1e-6)

})

