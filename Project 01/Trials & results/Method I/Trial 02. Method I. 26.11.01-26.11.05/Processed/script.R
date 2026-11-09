a <- read.csv("20261105-Amp-Lane1-22-Raw signals.csv")$Area /
  read.csv("20261105-Amp-LaneM-Raw signals.csv")$Area[1]
b <- read.csv("20261105-Kan-Lane1-22-Raw signals.csv")$Area /
  read.csv("20261105-Kan-LaneM-Raw signals.csv")$Area[1]

write.csv(
  data.frame(
    Lane = seq_along(a),
    "Amp normalized signal" = a,
    "Kan normalized signal" = b,
    check.names = FALSE
  ),
  "20261105-Normalized lane signals.csv",
  row.names = FALSE
)