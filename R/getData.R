#' extract data used in the figure
#'
#' \code{getData} extracts data used to plot different elements in the figure,
#' such as the heatmap, the row (column) names of the heatmap, etc.
#'
#' @param tree_hm the output of \code{\link{TreeHeatmap}}.
#' @param type It's chosen from \code{"heatmap", "row_name", "column_name",
#'   "title", "column_anno", "column_order"}.
#'
#' @export
#' @author Ruizhu Huang
#' @examples
#' library(TreeSummarizedExperiment)
#' library(ggtree)
#' library(ggplot2)
#' library(ggnewscale)
#' library(viridis)
#' library(dplyr)
#'
#' data(tinyTree)
#'
#' p1 <- c(rep(0.1/3, 3), rep(0.4/2, 2), rep(0.1, 5))
#' p2 <- c(rep(0.4/3, 3), rep(0.1/2, 2), rep(0.1, 5))
#' set.seed(1)
#' ct0 <- cbind(rmultinom(n = 5, size = 50, prob = p1),
#'              rmultinom(n = 5, size =50, prob = p2))
#' colnames(ct0) <- paste("S", 1:10, sep = "")
#' rownames(ct0) <- transNode(tree = tinyTree, node = 1:10)
#' oo <- sample(1:10)
#' ct0 <- ct0[, oo]
#'
#' ct <- rbind(colSums(ct0[1:3, ]),
#'             colSums(ct0[4:5, ]),
#'             ct0[6:10, ])
#' colnames(ct) <- paste("S", 1:10, sep = "")
#' rownames(ct) <- transNode(tree = tinyTree, node = c(13, 18, 6:10))
#'
#'
#' col_split <- ifelse(colnames(ct) %in% paste0("S", 1:5),
#'                     "A", "B")
#' names(col_split) <- colnames(ct)
#' # prepare the tree figure
#' tree_fig <- ggtree(tinyTree,
#'                    branch.length = "none",
#'                    layout = "rectangular", open.angle = 100) +
#'     #geom_text2(aes(label = label)) +
#'     geom_hilight(node = 18, fill = "orange", alpha = 0.3) +
#'     geom_hilight(node = 13, fill = "blue", alpha = 0.3)
#'  fig <- TreeHeatmap(tree = tinyTree, tree_fig = tree_fig,
#'             hm_data = ct, cluster_column = TRUE,
#'             column_split = col_split,
#'             column_anno = col_split,
#'             column_anno_gap = 0.6,
#'             column_anno_color = c("A" = "red", "B"= "blue"),
#'             show_colnames = TRUE,
#'             colnames_position = "bottom",
#'             colnames_angle = 90, colnames_size = 2,
#'             colnames_offset_y = -0.2,
#'             show_title = TRUE,
#'             title_offset_y = 1.5,
#'             title_color = "blue")
#' fig
#'
#' df_hm <- getData(tree_hm = fig, type = "heatmap")
#'
#' ct <- df_hm %>%
#'     dplyr::select(x, width, variable) %>%
#'     distinct()
#'
#' # generate data to do column annotation
#' set.seed(1)
#' ann <- matrix(sample(LETTERS[1:2], size = 3 * ncol(df_hm),
#'                      replace = TRUE),
#'               nrow = 3)
#' rownames(ann) <- paste0("g", 1:3)
#' colnames(ann) <- ct$variable
#' ann <- data.frame(ann) %>%
#'        mutate(y = min(df_hm$y) - seq_len(nrow(ann)),
#'        label = rownames(ann))
#' df_ann <- tidyr::gather(ann, variable, value, -c(y, label)) %>%
#'     left_join(ct)
#'
#' fig +
#'     new_scale_fill() +
#'     geom_tile(data = df_ann, aes(x, y-0.5,
#'                                  width = width, fill = value)) +
#'     scale_fill_viridis_d() +
#'     geom_text(data = df_ann, aes(x = min(x) - 1, y = y - 0.5,
#'                 label = label))


getData <- function(tree_hm, type = c("heatmap", "row_name", "column_name",
                                      "title", "column_anno", "column_order",
                                      "column_split")) {
    type <- match.arg(type)
    if (type == "heatmap") {
        out <- tree_hm$temp_data$hm_data
    }
    if (type == "row_name") {
        out <- tree_hm$temp_data$row_name
    }

    if (type == "column_name") {
        out <- tree_hm$temp_data$column_name
    }

    if (type == "title") {
        out <- tree_hm$temp_data$hm_title
    }

    if (type == "column_anno") {
        out <- tree_hm$temp_data$column_anno
    }

    if (type == "column_order") {
        out <- tree_hm$temp_data$column_order
    }

    if (type == "column_split") {
        out <- tree_hm$temp_data$column_split
    }
    return(out)
}
