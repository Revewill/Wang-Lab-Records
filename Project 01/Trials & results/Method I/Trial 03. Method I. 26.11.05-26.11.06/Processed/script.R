a <- read.csv("20261106-Amp-Lane1-13-Raw signals.csv")$Area /
  read.csv("20261106-Amp-LaneM-Raw signals.csv")$Area[1]
b <- read.csv("20261106-Kan-Lane1-13-Raw signals.csv")$Area /
  read.csv("20261106-Kan-LaneM-Raw signals.csv")$Area[1]

write.csv(
  data.frame(
    Lane = seq_along(a),
    "AmpR normalized signal" = a,
    "KanR normalized signal" = b,
    check.names = FALSE
  ),
  "20261106-Normalized lane signals.csv",
  row.names = FALSE
)