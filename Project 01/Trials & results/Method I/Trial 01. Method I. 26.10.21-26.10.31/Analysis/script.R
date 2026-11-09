# 首次运行时安装：
# install.packages("ggplot2")
library(ggplot2)

# ---------- 读取文件 ----------
message("请选择 NT 文件：result (2).txt")
nt_file <- file.choose()

message("请选择 T 文件：result.txt")
t_file <- file.choose()

read_result <- function(path) {
  d <- read.delim(
    path,
    header = TRUE,
    sep = "\t",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  names(d) <- sub("^#", "", names(d))
  required <- c(
    "Bulge type", "crRNA", "DNA", "Chromosome",
    "Position", "Direction", "Mismatches", "Bulge Size"
  )

  if (!all(required %in% names(d))) {
    stop("文件列名不符合 Cas-OFFinder 输出格式。")
  }

  # 无 bulge 排在前面；其后按坐标、错配数排列。
  # 保留同一坐标的不同对齐方式。
  d <- d[order(
    d[["Bulge Size"]],
    d[["Position"]],
    d[["Mismatches"]]
  ), ]

  rownames(d) <- NULL
  d
}

nt <- read_result(nt_file)
tt <- read_result(t_file)

# ---------- 作图函数 ----------
plot_alignment <- function(d, group) {

  n <- nrow(d)
  if (n == 0) stop("文件中没有命中结果。")

  max_width <- max(nchar(d$crRNA), nchar(d$DNA))
  cells <- vector("list", n)
  annotations <- vector("list", n)

  for (i in seq_len(n)) {

    query <- strsplit(toupper(d$crRNA[i]), "")[[1]]
    hit   <- strsplit(toupper(d$DNA[i]), "")[[1]]

    if (length(query) != length(hit)) {
      stop("第 ", i, " 行两条对齐序列长度不一致。")
    }

    len <- length(query)

    # 右对齐，使所有候选位点的 PAM 位于相同位置。
    x <- seq_len(len) + max_width - len
    pam <- seq_len(len) > len - 3

    dna_bulge <- query == "-" & !pam
    rna_bulge <- hit == "-" & !pam
    mismatch <- query != hit &
      query != "-" & hit != "-" & !pam

    state <- rep("Match", len)
    state[mismatch] <- "Mismatch"
    state[dna_bulge | rna_bulge] <- "Bulge"
    state[pam] <- "PAM"

    # 各位点之间留出空隙。
    y_top <- (n - i) * 3.15 + 1.05
    y_bottom <- y_top - 0.95

    cells[[i]] <- rbind(
      data.frame(
        x = x, y = y_top,
        base = query,
        state = ifelse(
          pam, "PAM",
          ifelse(dna_bulge | rna_bulge, "Bulge", "Match")
        ),
        row = "Query"
      ),
      data.frame(
        x = x, y = y_bottom,
        base = hit,
        state = state,
        row = "Genome"
      )
    )

    # 仅作为图内编号，同一坐标的对齐共享同一位点编号。
    locus_key <- paste(d$Chromosome, d$Position, d$Direction)
    locus_id <- match(locus_key[i], unique(locus_key))

    bulge_label <- if (d[["Bulge Size"]][i] == 0) {
      "No bulge"
    } else {
      paste0(
        d[["Bulge Size"]][i], "-nt ",
        d[["Bulge type"]][i], " bulge"
      )
    }

    annotations[[i]] <- data.frame(
      y_top = y_top,
      y_bottom = y_bottom,
      y_mid = (y_top + y_bottom) / 2,
      site = sprintf("%s-%02d", group, locus_id),
      coordinate = paste0(
        format(d$Position[i], big.mark = ",", scientific = FALSE),
        " (", d$Direction[i], ")"
      ),
      details = paste0(
        d$Mismatches[i], " mismatches",
        "\n", bulge_label
      )
    )
  }

  cells <- do.call(rbind, cells)
  ann <- do.call(rbind, annotations)

  # PAM 前增加少量空白。
  cells$x[cells$x > max_width - 3] <-
    cells$x[cells$x > max_width - 3] + 0.45

  fills <- c(
    Match = "#F5F6F7",
    Mismatch = "#F5D6D3",
    Bulge = "#D7E8F5",
    PAM = "#E4E7EB"
  )

  text_colors <- c(
    Match = "#35404A",
    Mismatch = "#B52E29",
    Bulge = "#216492",
    PAM = "#35404A"
  )

  cells$state <- factor(
    cells$state,
    levels = c("Match", "Mismatch", "Bulge", "PAM")
  )

  group_color <- if (group == "NT") "#245B83" else "#B76560"

  p <- ggplot(cells, aes(x, y)) +

    geom_tile(
      aes(fill = state),
      width = 0.91,
      height = 0.84,
      show.legend = TRUE
    ) +

    geom_text(
      aes(label = base, color = state),
      family = "mono",
      size = 4.0,
      show.legend = FALSE
    ) +

    # 左侧位点编号及 Query / Genome 行标签。
    geom_text(
      data = ann,
      aes(x = -6.4, y = y_mid, label = site),
      inherit.aes = FALSE,
      hjust = 0,
      fontface = "bold",
      color = group_color,
      size = 3.8
    ) +

    geom_text(
      data = ann,
      aes(x = -1.0, y = y_top),
      label = "Query",
      inherit.aes = FALSE,
      hjust = 1,
      color = "#46515A",
      size = 3.3
    ) +

    geom_text(
      data = ann,
      aes(x = -1.0, y = y_bottom),
      label = "Genome",
      inherit.aes = FALSE,
      hjust = 1,
      color = "#46515A",
      size = 3.3
    ) +

    # 右侧基因组坐标及错配 / bulge 信息。
    geom_text(
      data = ann,
      aes(
        x = max_width + 2.0,
        y = y_top,
        label = coordinate
      ),
      inherit.aes = FALSE,
      hjust = 0,
      size = 3.4,
      color = "#35404A"
    ) +

    geom_text(
      data = ann,
      aes(
        x = max_width + 2.0,
        y = y_bottom - 0.12,
        label = details
      ),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 0.5,
      lineheight = 0.95,
      size = 3.0,
      color = "#59636B"
    ) +

    scale_fill_manual(
      values = fills,
      breaks = c("Mismatch", "Bulge", "PAM"),
      drop = FALSE,
      name = NULL
    ) +

    scale_color_manual(values = text_colors, guide = "none") +

    scale_x_continuous(
      limits = c(-6.7, max_width + 10),
      expand = c(0, 0)
    ) +

    scale_y_continuous(
      limits = c(-1.2, max(ann$y_top) + 1),
      expand = c(0, 0)
    ) +

    labs(title = paste0(group, " | Candidate sequence alignments")) +

    theme_void(base_family = "sans") +

    theme(
      plot.title = element_text(
        size = 17,
        face = "bold",
        color = group_color,
        margin = margin(b = 15)
      ),
      legend.position = "top",
      legend.justification = "left",
      legend.text = element_text(size = 10),
      legend.key.size = grid::unit(0.38, "cm"),
      legend.margin = margin(b = 8),
      plot.margin = margin(15, 20, 15, 15),
      plot.background = element_rect(
        fill = "white",
        color = NA
      )
    )

  p
}

# ---------- 显示并导出 ----------
p_nt <- plot_alignment(nt, "NT")
p_t  <- plot_alignment(tt, "T")

print(p_nt)
print(p_t)

out_dir <- file.path(dirname(nt_file), "Cas_OFFinder_alignment")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

save_alignment <- function(p, d, group) {

  # 根据命中条数调整高度，避免文字挤在一起。
  height <- 1.5 + 0.66 * nrow(d)

  ggsave(
    filename = file.path(out_dir, paste0(group, "_alignment.png")),
    plot = p,
    width = 13,
    height = height,
    units = "in",
    dpi = 400,
    bg = "white",
    limitsize = FALSE
  )

  ggsave(
    filename = file.path(out_dir, paste0(group, "_alignment.pdf")),
    plot = p,
    width = 13,
    height = height,
    units = "in",
    bg = "white",
    limitsize = FALSE
  )
}

save_alignment(p_nt, nt, "NT")
save_alignment(p_t, tt, "T")

message("输出位置：", normalizePath(out_dir))