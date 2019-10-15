#' Heatmap at arbitary levels of a tree
#'
#' \code{TreeHeatmap} displays a heatmap at a arbitary level of a tree.
#'
#' @param tree A phylo object
#' @param tree_fig A ggtree object that outputs from
#'   \code{\link[ggtree]{ggtree}}.
#' @param hm_data A data frame. Data to plot heatmap. Its rownames should be
#'   able to match to nodes of the \strong{tree}.
#' @param tree_hm_gap A numeric value to specify the gap between the tree and
#'   the heatmap.
#' @param rel_width A numeric value to specify the width of heatmap relative to
#'   the width of the tree. For example, \code{rel_width = 1}, the width of
#'   heatmap is the same as the width of the tree.
#' @param column_order A character vector that includes the column name of
#'   \strong{hm_data} to specify the display order of the heatmap. It's ignored
#'   when \strong{column_split} is provided.
#' @param column_split A named character vector that gives the group information
#'   about columns to split the heatmap. It's named by the colnames of
#'   \strong{hm_data}.
#' @param column_split_gap A numeric value to specify the gap between the
#'   columns of heatmap that are split.
#' @param column_anno A named vector to specify labels that are used to
#'   annotate columns of heatmap.
#' @param column_anno_size A numeric value to specify the size of the annotation
#'   bar
#' @param column_anno_color A named vector to specify colors that are used to
#'   annotate columns of heatmap.
#' @param column_anno_gap A numeric value to specify the gap between the
#'   column annotation bar and the heatmap.
#' @param legend_title_hm The legend title of the heatmap.
#' @param legend_title_column_anno The legend title of the column annotation.
#' @param show_colnames A logical value to specify whether column names should
#'   be displayed. The default is FALSE.
#' @param colnames_position The position of column names, "top" or "bottom".
#' @param colnames_angle A numeric value. The angle of column names.
#' @param colnames_offset_x A numeric value to shift column names on x-axis. The
#'   defaut is 0.
#' @param colnames_offset_y A numeric value to shift column names on y-axis. The
#'   defaut is 0.
#' @param colnames_size A numeric value to specify the size of column names.
#' @param colnames_hjust The hjust for column names: 0 (left aligned); 0.5
#'   (centered); 1 (right aligned).
#' @param show_rownames A logical value to specify whether row names should
#'   be displayed. The default is FALSE.
#' @param rownames_angle A numeric value. The angle of row names.
#' @param rownames_offset_x A numeric value to shift row names on x-axis. The
#'   defaut is 0.
#' @param rownames_offset_y A numeric value to shift row names on y-axis. The
#'   defaut is 0.
#' @param rownames_size A numeric value to specify the size of row names.
#' @param rownames_hjust The hjust for row names: 0 (left aligned); 0.5
#'   (centered); 1 (right aligned).
#' @param show_title A logical value to specify whether the title should
#'   be displayed. The default is FALSE.
#' @param title_hm The title of heatmap
#' @param title_fontface The fontface of the title. The default is "bold".
#' @param title_color  The color of the title. The default is "black".
#' @param title_size The size of the title. The default is 3.
#' @param title_angle The angle of the title. The default is 0.
#' @param title_offset_x A numeric value to shift title along x-axis. The
#'   defaut is 0.
#' @param title_offset_y A numeric value to shift title along y-axis. The
#'   defaut is 2.
#' @param title_hjust The hjust for title: 0 (left aligned); 0.5
#'   (centered); 1 (right aligned). The default is 0.5
#' @param cluster_column A logical value, TRUE or FALSE. It specifies whether
#'   columns of the heatmap should be ordered by similarity. The default is
#'   TRUE. This is ignored when \strong{column_order} is given.
#' @importFrom TreeSummarizedExperiment transNode findOS
#' @importFrom ggtree ggtree
#' @importFrom tidyr gather
#' @importFrom dplyr mutate select distinct "%>%"
#' @importFrom ggplot2 geom_tile geom_segment scale_color_manual labs geom_text scale_fill_viridis_c aes
#' @importFrom ggnewscale new_scale_color
#' @importFrom viridis viridis
#' @importFrom stats hclust dist
#' @export
#' @author Ruizhu Huang
#' @examples
#'
#' library(TreeSummarizedExperiment)
#' library(ggtree)
#' library(ggplot2)
#' library(scales)
#'
#' data(tinyTree)
#' p1 <- c(0.1, 0.4, rep(0.1, 5))
#' p2 <- c(0.4, 0.1, rep(0.1, 5))
#' set.seed(1)
#' ct <- cbind(rmultinom(n = 5, size = 50, prob = p1),
#'             rmultinom(n = 5, size =50, prob = p2))
#' colnames(ct) <- paste("S", 1:10, sep = "")
#' rownames(ct) <- transNode(tree = tinyTree, node = c(13, 18, 6:10))
#'
#' oo <- sample(1:10)
#' ct <- ct[, oo]
#' # prepare the tree figure
#' tree_fig <- ggtree(tinyTree,
#'                    branch.length = "none",
#'                    layout = "rectangular") +
#'     #geom_text2(aes(label = label)) +
#'     geom_hilight(node = 18, fill = "orange", alpha = 0.3) +
#'     geom_hilight(node = 13, fill = "blue", alpha = 0.3)
#' # figure 1
#' TreeHeatmap(tree = tinyTree, tree_fig = tree_fig, hm_data = ct)
#'
#' # figure 2: order column by similarity
#' TreeHeatmap(tree = tinyTree, tree_fig = tree_fig,
#'             hm_data = ct, cluster_column = TRUE)
#'
#'  # figure 3: split columns
#' col_split <- ifelse(colnames(ct) %in% paste0("S", 1:5),
#'                     "A", "B")
#' names(col_split) <- colnames(ct)
#' TreeHeatmap(tree = tinyTree, tree_fig = tree_fig,
#'             hm_data = ct, cluster_column = TRUE,
#'             column_split = col_split)
#' # figure 4: annotate columns
#' col_anno <- col_split
#' TreeHeatmap(tree = tinyTree, tree_fig = tree_fig,
#'             hm_data = ct, cluster_column = TRUE,
#'             column_split = col_split,
#'             column_anno = col_anno,
#'             column_anno_gap = 0.5)
#' # figure 5: change annotation colors
#' TreeHeatmap(tree = tinyTree, tree_fig = tree_fig,
#'             hm_data = ct, cluster_column = TRUE,
#'             column_split = col_split,
#'             column_anno = col_anno,
#'             column_anno_gap = 0.6,
#'             column_anno_color = c("A" = "red", "B"= "blue"))
#'
#' # figure 6: add colnames
#' TreeHeatmap(tree = tinyTree, tree_fig = tree_fig,
#'             hm_data = ct, cluster_column = TRUE,
#'             column_split = col_split,
#'             column_anno = col_anno,
#'             column_anno_gap = 0.6,
#'             column_anno_color = c("A" = "red", "B"= "blue"),
#'             show_colnames = TRUE,
#'             colnames_position = "bottom",
#'             colnames_angle = 90, colnames_size = 2,
#'             colnames_offset_y = -0.2)
#'
#' # figure 7: add title
#'  fig <- TreeHeatmap(tree = tinyTree, tree_fig = tree_fig,
#'             hm_data = ct, cluster_column = TRUE,
#'             column_split = col_split,
#'             column_anno = col_anno,
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
#' # use expand_limits to display the missing part of row names
#' (fig <- fig + expand_limits(x = c(0, 15)))
#'
#' # change colors
#' fig +
#' scale_fill_viridis_c(option = "D")
#'
#' fig +
#'  scale_fill_gradientn(colours = c("blue","yellow","red"),
#'                       values = scales::rescale(c(5, 8, 10)),
#'                       guide = "colorbar", limits=c(5, 10))

TreeHeatmap <- function(tree, tree_fig, hm_data,
                        tree_hm_gap = 0,
                        rel_width = 1,
                        column_order = NULL,
                        column_split = NULL,
                        column_anno = NULL,
                        column_anno_size = 1,
                        column_split_gap = 0.2,
                        column_anno_color = NULL,
                        column_anno_gap = 0.1,
                        legend_title_hm = "Expression",
                        legend_title_column_anno = "group",
                        show_colnames = FALSE,
                        colnames_position = "top",
                        colnames_angle = 0,
                        colnames_offset_x = 0,
                        colnames_offset_y = 0,
                        colnames_size = 4,
                        colnames_hjust = 0.5,
                        show_rownames = FALSE,
                        rownames_angle = 0,
                        rownames_offset_x = 0,
                        rownames_offset_y = 0,
                        rownames_size = 4,
                        rownames_hjust = 0.5,
                        show_title = FALSE,
                        title_hm = "First heatmap",
                        title_fontface = "bold",
                        title_color = "black",
                        title_size = 3,
                        title_angle = 0,
                        title_offset_x = 0,
                        title_offset_y = 2,
                        title_hjust = 0.5,
                        cluster_column = FALSE){


    if (!is.null(column_order)) {
        if (!all(column_order %in% colnames(hm_data))) {
            stop("column_order: Some columns don't exist in hm_data.")
        }
    }

    if (!is.null(column_split)) {
        if (is.null(names(column_split))) {
            stop("column_split: should be named by colnames of hm_data")
        }
        if (!all(names(column_split) %in% colnames(hm_data))) {
            stop("column_split: Some columns don't exist in hm_data.")
        }
    }

    if (!is.null(column_anno_color)) {
        if (is.null(names(column_anno_color))) {
            stop("column_anno_color: should be named by column_anno")
        }

    }

    ## tree: data

    df <- tree_fig$data

    # ------------------------ heatmap -----------------------
    # data
    hm_df <- data.frame(hm_data, check.names = FALSE)

    # heatmap: node
    rnam_hm <- rownames(hm_df)
    node_hm <- transNode(tree = tree, node = rnam_hm)
    hm_df$node <- node_hm
    hm_df$row_label <- rnam_hm

    # heatmap: y
    desd_hm <- findOS(tree = tree, node = node_hm,
                      only.leaf = FALSE, self.include = TRUE)
    y_hm <- lapply(desd_hm, FUN = function(x){
        xx <- match(x, df$node)
        y <- df$y[xx]
        # the middle point
        mean(range(y, na.rm = TRUE))
    })
    hm_df$y <- unlist(y_hm)

    # heatmap: height of a row
    h_hm <- lapply(desd_hm, FUN = function(x){
        xx <- match(x, df$node)
        y <- df$y[xx]

        cy <- colnames(df)
        if ("scale" %in% cy) {
            dt <- unique(df$scale[xx])
        } else {
            dt <- 1
        }
        # the distance
        diff(range(y, na.rm = TRUE)) + dt
    })
    hm_df$height <- unlist(h_hm)

    # heatmap: width of a column
    width_hm <- rel_width * (df$x %>%
                                 range(na.rm = TRUE) %>%
                                 diff)/ncol(hm_data)
    hm_df$width <- width_hm

    # heatmap: long form
    variable <- node <- row_label <- NULL
    value <- x <- y <- height <- width <- NULL
    hm_dt <- gather(hm_df, variable, value,
                    -c(y, node, row_label,
                       width, height))

    # heatmap: column order
    if (!is.null(column_split)) {
        # 1. column_split is given, ignore column order. The order within the
        #    same slice is determined by the input order of column split
        # 2. column_split is given and column_cluster = TRUE, the order within
        #    the same slice is determined by the column similarity
        column_split <- factor(column_split, levels = unique(column_split))
        if (cluster_column) {
            # column similarity within the same slice
            ind_split <- lapply(levels(column_split), FUN = function(x){
                names(column_split)[column_split == x]
            })
            similar <- lapply(ind_split, FUN = function(x) {
                if (length(x) > 2) {
                    st <- hm_data[, x, drop = FALSE]
                    xx <- hclust(dist(t(st)), "ave")
                    colnames(st)[xx$order]
                } else { x }
            })
            similar_order <- unlist(similar)
            column_split <- column_split[similar_order]
        }
        split_level <- sort(column_split)

        if (!is.null(column_order)) {
            warnings("column_order is ignored when column_split is given")}
        column_order <- names(split_level)
    } else {
        ## column_split isn't given
        # 1. column_order is given, use column_order and ignore column
        #    similarity.
        # 2. column_order isn't given but allow cluster columns, order columns
        #    by similarity
        # 3. column_order isn't given and cluter columns is not allowed, use the
        #    original column order.

        split_level <- rep(0, ncol(hm_data))
        names(split_level) <- colnames(hm_data)
        split_level <- factor(split_level)
        if (is.null(column_order) & !cluster_column) {
            column_order <- colnames(hm_data)}
        if (is.null(column_order) & cluster_column) {
            hc <- hclust(dist(t(hm_data)), "ave")
            column_order <- colnames(hm_data)[hc$order]}
        if (!is.null(column_order) & cluster_column) {
            warnings("cluster_column is ignored because column_order is given")
        }
    }

    # heatmap: x
    hm_dt <- hm_dt %>%
        mutate(variable = factor(variable, levels = column_order)) %>%
        mutate(column_order = as.numeric(variable) -1) %>%
        mutate(split_level = split_level[variable]) %>%
        mutate(split_level = as.numeric(split_level) - 1) %>%
        mutate(x = max(df$x, na.rm = TRUE) +
                   tree_hm_gap + width_hm/2 +
                   column_order  * width_hm +
                   split_level * column_split_gap) %>%
        select(node, row_label, x, y,
               height, width, variable,
               value, column_order, split_level)

    # # -------------------- tree + heatmap --------------------
    p <- tree_fig +
        geom_tile(data = hm_dt,
                  aes(x = x, y = y,
                      height = height,
                      fill = value,
                      width = width),
                  inherit.aes = FALSE)
    p <- p + scale_fill_viridis_c()
    p

    # # -------------------- heatmap annotation ----------------
    if (!is.null(column_anno)) {
        if (is.null(column_anno_color)) {
            anno_uc <- unique(column_anno)
            column_anno_color <- viridis(length(anno_uc))
            names(column_anno_color) <- anno_uc
        }
        # column annotations
        xend <- yend <- anno_group <- anno_color <- NULL
        anno_df <- hm_dt %>%
            select(variable, x, width) %>%
            distinct() %>%
            mutate(
                variable = as.character(variable),
                x = x - 0.5 * width,
                xend = x + width,
                y = max(df$y, na.rm = TRUE) +
                    column_anno_gap,
                yend = max(df$y, na.rm = TRUE) +
                    column_anno_gap,
                anno_group = column_anno[variable],
                anno_color = column_anno_color[anno_group])
        anno_color <- anno_df$anno_color
        names(anno_color) <- anno_df$anno_group


        p <- p +
            new_scale_color()+
            geom_segment(data = anno_df,
                         aes(x = x, y = y,
                             xend = xend,
                             yend = yend,
                             color = anno_group),
                         inherit.aes = FALSE,
                         size = column_anno_size) +
            scale_color_manual(values = anno_color) +
            labs(fill = legend_title_hm,
                 color = legend_title_column_anno)
    }

    # # -------------------- heatmap column names ----------------
    if (show_colnames) {
        y_top <- y_bottom <- NULL
        cn_df <- hm_dt %>%
            select(variable, x, width) %>%
            distinct() %>%
            mutate(
                variable = as.character(variable),
                y_top = max(hm_dt$y + 0.5 * hm_dt$height) +
                    column_anno_gap,
                y_bottom = min(hm_dt$y - 0.5* hm_dt$height),
                y = ifelse(colnames_position == "top",
                           y_top, y_bottom))
        p <- p + geom_text(data = cn_df,
                           aes(x = x, y = y,
                               label = variable),
                           size = colnames_size, inherit.aes = FALSE,
                           angle = colnames_angle,
                           nudge_x = colnames_offset_x,
                           nudge_y = colnames_offset_y,
                           hjust = colnames_hjust)
    }

    if (show_rownames) {
        rn_df <- hm_dt %>%
            select(y, width, row_label) %>%
            distinct() %>%
            mutate(x = max(hm_dt$x + 0.5*hm_dt$width))
        p <- p + geom_text(data = rn_df,
                           aes(x = x, y = y,
                               label = row_label),
                           size = rownames_size, inherit.aes = FALSE,
                           angle = rownames_angle,
                           nudge_x = rownames_offset_x,
                           nudge_y = rownames_offset_y,
                           hjust = rownames_hjust)
    }

    # # -------------------- heatmap title ----------------
    if (show_title) {
        label <- NULL
        title_df <- data.frame(x = mean(range(hm_dt$x, na.rm = TRUE)),
                               y = max(hm_dt$y),
                               label = title_hm)
        p <- p + geom_text(data = title_df,
                           aes(x = x, y = y,
                               label = label),
                           inherit.aes = FALSE,
                           fontface = title_fontface,
                           colour = title_color,
                           size = title_size,
                           angle = title_angle,
                           nudge_x = title_offset_x,
                           nudge_y = title_offset_y,
                           hjust = title_hjust)
    }
    return(p)
}
