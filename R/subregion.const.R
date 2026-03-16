#' Domain decomposition based on a triangulated mesh
#'
#' @description
#' `subregion_const()` constructs a **non-overlapping domain decomposition**
#' of a triangulated spatial domain. The function partitions the triangulation
#' into disjoint sub-triangulations, each representing a spatial subregion.
#'
#' The decomposition is based on the connectivity structure of the triangulation.
#' Around each triangle, the algorithm constructs a local neighborhood using
#' adjacency layers, and then selects a subset of non-overlapping neighborhoods
#' as seeds for the final partition. Remaining triangles are subsequently assigned
#' to the seed whose neighborhood overlaps the most.
#'
#' Let \eqn{\triangle} denote the triangulation of the spatial domain
#' \eqn{\Omega}. The algorithm partitions \eqn{\triangle} into
#' \eqn{m} disjoint sub-triangulations
#' \eqn{\triangle_1, \ldots, \triangle_m} such that
#'
#' \deqn{\triangle = \cup_{k=1}^m \triangle_k,}
#'
#' with \eqn{\triangle_k \cap \triangle_{k'} = \emptyset} for \eqn{k \ne k'}.
#' Each sub-triangulation \eqn{\triangle_k} corresponds to a spatial subregion
#' \eqn{\omega_k}, providing a mesh-based representation of the domain partition.
#'
#' This decomposition is used in distributed inference procedures such as
#' **DISQ**, where each subregion can be processed on a separate computing node.
#'
#' @param Tr Integer matrix with 3 columns. Each row contains the vertex
#'   indices defining one triangle in the triangulation.
#' @param V Numeric matrix with 2 columns giving the vertex coordinates
#'   \eqn{(x, y)}.
#' @param nt Integer. Total number of triangles. Default is `nrow(Tr)`.
#' @param n.layer Positive integer controlling the neighborhood size used
#'   when constructing candidate subregions. Larger values produce larger
#'   candidate neighborhoods and therefore coarser partitions.
#'
#' @details
#' The partitioning algorithm proceeds as follows:
#'
#' **Step 1. Construct candidate neighborhoods.**
#'
#' For each triangle \eqn{l}, construct a neighborhood subregion
#' \eqn{R_l} consisting of triangles reachable within `n.layer`
#' adjacency layers in the triangulation graph.
#'
#' **Step 2. Select disjoint seed subregions.**
#'
#' From the candidate neighborhoods \eqn{R_l}, select a subset of seeds
#' whose neighborhoods are mutually disjoint. These seeds form the
#' initial collection of subregions.
#'
#' **Step 3. Assign remaining triangles.**
#'
#' For each triangle not yet assigned to any subregion, compute its
#' overlap with the existing seed neighborhoods and assign it to the
#' subregion with the largest overlap.
#'
#' The final output is a set of disjoint triangle index sets representing
#' the sub-triangulations.
#'
#' The implementation assumes the following helper functions are available:
#'
#' - `tdata(V, Tr)$TV` returns triangle-vertex adjacency information.
#' - `ring(Tr0, V, Tr, TV, n.layer)` constructs a neighborhood around a
#'   triangle using adjacency layers and returns the triangles contained
#'   in that neighborhood.
#'
#' @return
#' A list `triangles` where each element is an integer vector containing
#' the row indices of `Tr` corresponding to triangles in one subregion.
#'
#' @examples
#' \dontrun{
#' ## Step 1: Create a triangulated mesh
#' boundaries <- matrix(
#'   c(0, 0,
#'     0, 1,
#'     1, 1,
#'     1, 0),
#'   nrow = 4, byrow = TRUE
#' )
#'
#' tri <- TriMesh(boundaries, n = 13)
#' Tr  <- tri$Tr
#' V   <- round(tri$V, 3)
#'
#' ## Step 2: Perform domain decomposition
#' triangles <- subregion_const(Tr = Tr, V = V, n.layer = 1)
#'
#' ## Step 3: Visualize the resulting subregions
#' sub_tri_plot(triangles, V = V, Tr = Tr, main = expression(Delta))
#' }
#'
#' @export
subregion_const <- function(nt = nrow(Tr), Tr, V, n.layer = 1) {

  stopifnot(is.matrix(Tr), ncol(Tr) == 3)
  stopifnot(is.matrix(V),  ncol(V)  == 2)
  stopifnot(is.numeric(n.layer), length(n.layer) == 1, n.layer >= 1)

  nt <- as.integer(nt)

  # --------------------------------------------
  # Step 1: Build candidate subregions
  # --------------------------------------------
  TV <- tdata(V, Tr)$TV

  sub_regions <- vector("list", length = nrow(Tr))
  for (j in seq_len(nrow(Tr))) {
    sub_regions[[j]] <- ring(
      Tr0 = Tr[j, ],
      V = V,
      Tr = Tr,
      TV = TV,
      n.layer = n.layer
    )
  }

  # Utility: map a set of triangles to row indices in Tr
  match_triangles <- function(T_sub, Tr_ref) {
    match(data.frame(t(T_sub)), data.frame(t(Tr_ref)))
  }

  # --------------------------------------------
  # Step 2: Select disjoint seed subregions
  # --------------------------------------------
  triangles <- vector("list", 0)
  collection <- integer(0)

  seed1_idx <- match_triangles(sub_regions[[1]]$T.sub, Tr)
  triangles[[1]] <- seed1_idx
  collection <- 1L

  for (ll in seq_len(nt)) {

    before <- length(triangles)

    pool <- sort(unique(unlist(triangles)))

    chosen_seed <- NA_integer_

    candidates <- setdiff(seq_len(nrow(Tr)), pool)
    if (length(candidates) == 0) break

    for (j in candidates) {

      tr_j <- match_triangles(sub_regions[[j]]$T.sub, Tr)
      if (anyNA(tr_j)) next

      ok_with_all <- TRUE
      for (k in seq_along(collection)) {
        tr_k <- match_triangles(sub_regions[[collection[k]]]$T.sub, Tr)
        if (anyNA(tr_k)) next
        if (length(intersect(tr_k, tr_j)) > 0) {
          ok_with_all <- FALSE
          break
        }
      }

      if (ok_with_all) {
        chosen_seed <- j
        break
      }
    }

    if (is.na(chosen_seed)) break

    triangles[[length(triangles) + 1L]] <-
      match_triangles(sub_regions[[chosen_seed]]$T.sub, Tr)

    collection <- c(collection, chosen_seed)

    after <- length(triangles)
    if (before == after) break
  }

  # --------------------------------------------
  # Step 3: Assign remaining triangles
  # --------------------------------------------
  pool <- sort(unique(unlist(triangles)))
  leftovers <- setdiff(seq_len(nrow(Tr)), pool)

  if (length(leftovers) > 0) {
    for (j in leftovers) {

      tr_j <- match_triangles(sub_regions[[j]]$T.sub, Tr)
      if (anyNA(tr_j)) next

      overlap_counts <- integer(length(triangles))
      for (k in seq_along(triangles)) {
        overlap_counts[k] <- length(intersect(triangles[[k]], tr_j))
      }

      kbest <- which.max(overlap_counts)
      triangles[[kbest]] <- sort(unique(c(triangles[[kbest]], j)))
    }
  }

  triangles
}


#' Plot a triangulation domain decomposition
#'
#' @description
#' `sub_tri_plot()` visualizes a triangulation partition produced by
#' `subregion_const()`. Each subregion is displayed using a different color,
#' allowing the user to inspect the resulting domain decomposition.
#'
#' Optionally, a designated **master subregion** can be highlighted
#' in red. By default, the largest subregion is selected.
#'
#' @param triangles List of triangle index sets. Each element corresponds
#'   to a subregion and contains row indices of `Tr`.
#' @param V Numeric matrix with 2 columns giving vertex coordinates.
#' @param Tr Integer matrix with 3 columns defining the triangulation
#'   connectivity.
#' @param main Plot title. Can be a character string or an expression.
#' @param colors Optional vector of colors used for the subregions.
#'   If `NULL`, colors are generated automatically.
#' @param master Integer index indicating the subregion to highlight.
#'   By default, the largest subregion is selected.
#' @param border Color used for triangle borders. Default is `"white"`.
#'
#' @return
#' Invisibly returns `NULL`.
#'
#' @examples
#' \dontrun{
#' boundaries <- matrix(
#'   c(0, 0,
#'     0, 1,
#'     1, 1,
#'     1, 0),
#'   nrow = 4, byrow = TRUE
#' )
#'
#' tri <- TriMesh(boundaries, n = 13)
#' Tr  <- tri$Tr
#' V   <- round(tri$V, 3)
#'
#' triangles <- subregion_const(Tr = Tr, V = V, n.layer = 1)
#'
#' sub_tri_plot(
#'   triangles = triangles,
#'   V = V,
#'   Tr = Tr,
#'   main = expression(Delta)
#' )
#' }
#'
#' @export
sub_tri_plot <- function(triangles, V, Tr,
                         main = expression(Delta),
                         colors = NULL,
                         master = which.max(vapply(triangles, length, numeric(1))),
                         border = "white") {

  stopifnot(is.list(triangles), length(triangles) >= 1)
  stopifnot(is.matrix(V), ncol(V) == 2)
  stopifnot(is.matrix(Tr), ncol(Tr) == 3)
  stopifnot(master >= 1, master <= length(triangles))

  xlim <- range(V[, 1])
  ylim <- range(V[, 2])

  plot(
    NA,
    xlim = xlim, ylim = ylim,
    xlab = "", ylab = "",
    axes = FALSE, asp = 1,
    main = main
  )

  if (is.null(colors)) {
    colors <- grDevices::hcl.colors(length(triangles), palette = "Set3")
  } else {
    stopifnot(length(colors) >= length(triangles))
  }

  for (i in seq_along(triangles)) {
    col_i <- if (i == master) "red" else colors[i]
    tri_ids <- triangles[[i]]

    for (j in seq_along(tri_ids)) {
      polygon(
        V[Tr[tri_ids[j], ], , drop = FALSE],
        border = border,
        col = col_i
      )
    }
  }

  invisible(NULL)
}
