context("Internal cluster validity functions")


test_that("Internal cluster validity measures work", {

    ft_mat <- data.frame(a = c(0, 0, 1, 1),
                         b = c(1, 1, 0, 0),
                         c = c(1, 1, 1, 1))

    d_mat <- dist(ft_mat)
    hc <- hclust(d_mat)
    clusters <- cutree(hc, 2)

    expect_equal(avgSilhouette(clusters, d_mat), 1)

    expect_equal(internalClusterValidity(2, hc, d_mat),
                 data.frame(k = 2,
                            Average_Silhouette = 1))

    expect_equal(clusterValidityPerK(ft_mat),
                 data.frame(k = c(2, 3),
                            Average_Silhouette = c(1.0, 0.5)))

})
