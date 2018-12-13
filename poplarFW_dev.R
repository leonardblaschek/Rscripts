library(dplyr)
library(tidyr)
library(broom)
library(purrr)
library(ggplot2)

poplar <-
  read.csv("file:///home/leonard/Documents/Uni/Phloroglucinol/poplar_foodweb.csv")
poplar <- poplar[, -16]

poplar$cell.type <- recode(poplar$cell.type,
  "F" = "Fibre",
  "V" = "Vessel",
  "R" = "Ray",
  "CB" = "Cambium"
)
poplar$adj.cell.type <- recode(poplar$adj.cell.type,
  "F" = "Fibre",
  "V" = "Vessel",
  "R" = "Ray",
  "CB" = "Cambium",
  "PA" = "Parenchyma"
)

poplar$Distance <-
  apply(poplar[, c("X", "Y", "ref.x1", "ref.x2", "ref.y1", "ref.y2")],
        1 ,
        function(x) {
          a <- c(x[1], x[2])
          b <- c((x[3] / 5.9), (x[5] / 5.9))
          c <- c((x[4] / 5.9), (x[6] / 5.9))
          v1 <- b - c
          v2 <- a - b
          m <- cbind(v1, v2)
          d <- abs(det(m)) / sqrt(sum(v1 * v1))
          d
        })

poplar$diff <- poplar$OD.stained - poplar$OD.unstained
poplar.bg <-
  subset(
    poplar,
    cell.type == "Cambium" & adj.cell.type == "Cambium",
    select = c("genotype", "replicate", "technical", "diff")
  )
colnames(poplar.bg)[4] <- "OD.bg"

poplar.bg <- poplar.bg %>%
  group_by(genotype, replicate, technical) %>%
  mutate(OD.bg = mean(OD.bg, na.rm = TRUE))

poplar <-
  full_join(poplar,
            unique(poplar.bg),
            by = c("genotype", "replicate", "technical"))
poplar$diff.adj <- poplar$diff - poplar$OD.bg

poplar.bin <- poplar %>%
  mutate(bin = cut(
    Distance,
    breaks = c(-Inf, 50, 100, Inf),
    labels = c("I", "II", "III")
  ))

poplar.bin.count <- poplar.bin %>%
  group_by(genotype, replicate, cell.type, adj.cell.type, bin) %>%
  summarise(count = n())

poplar.bin.pre <- poplar.bin %>%
  filter(cell.type %in% c("Vessel", "Fibre", "Ray"), adj.cell.type %in% c("Vessel", "Fibre", "Ray")) %>%
  group_by(genotype, bin, cell.type, adj.cell.type, replicate) %>%
  summarise(od = mean(diff.adj))

# poplar.bin.avg <- poplar.bin.pre %>%
#   group_by(genotype, bin, cell.type, adj.cell.type) %>%
#   summarise(od = mean(od))

poplar.bin.spread <- poplar.bin.pre %>%
  unite(cell.wall, cell.type, adj.cell.type) %>%
  spread(cell.wall, od)

poplar.lm <- ungroup(poplar.bin.spread) %>%
  select(-genotype, replicate) %>%
  group_by(bin) %>%
  do("V -> R" = lm(Ray_Vessel ~ Vessel_Vessel, data = .),
     "V -> F" = lm(Fibre_Vessel ~ Vessel_Vessel, data = .),
     "F -> V" = lm(Vessel_Fibre ~ Fibre_Fibre, data = .),
     "F -> R" = lm(Ray_Fibre ~ Fibre_Fibre, data = .),
     "R -> V" = lm(Vessel_Ray ~ Ray_Ray, data = .),
     "R -> F" = lm(Fibre_Ray ~ Ray_Ray, data = .)
     )
poplar.lm <- poplar.lm %>%
  gather(-bin, key ="relation", value = "regression") %>%
  unite(relation_bin, relation, bin)

poplar.lm.tidy <- map_dfr(poplar.lm$regression, glance, .id = "relation_bin")
poplar.lm.tidy <- full_join(
  map_dfr(poplar.lm$regression, glance, .id = "relation_bin"),
  subset(map_dfr(poplar.lm$regression, tidy, .id = "relation_bin"), term != "(Intercept)"),
  by = "relation_bin")
poplar.lm.tidy$relation_bin <- poplar.lm$relation_bin
write.csv(poplar.lm.tidy, file = "poplar_foodwebs.csv")

pdf("r_squared_genotypes.pdf", width = 15)
ggplot(data = poplar.lm.tidy, aes(x = relation_bin, y = adj.r.squared)) +
  geom_bar(stat = "identity")
dev.off()


