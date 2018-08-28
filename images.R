library(sysfonts)
library(showtext)
library(cowplot)
library(png)
library(grid)

# import Helvetica Neue
font_add("Helvetica",
         regular = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Lt.otf",
         italic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-LtIt.otf",
         bold = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-Bd.otf",
         bolditalic = "/prop_fonts/01. Helvetica     [1957 - Max Miedinger]/HelveticaNeueLTStd-BdIt.otf")
showtext_auto()

WTstained <- rasterGrob(readPNG("~/Documents/Uni/Phloroglucinol/18-06_draft/Images/stained/WT_stained.png"))
cl1stained <- rasterGrob(readPNG("~/Documents/Uni/Phloroglucinol/18-06_draft/Images/stained/4cl1_stained.png"))
cl2stained <- rasterGrob(readPNG("~/Documents/Uni/Phloroglucinol/18-06_draft/Images/stained/4cl2_stained.png"))
cl1x2stained <- rasterGrob(readPNG("~/Documents/Uni/Phloroglucinol/18-06_draft/Images/stained/4cl1x2_stained.png"))
ccoaomtstained <- rasterGrob(readPNG("~/Documents/Uni/Phloroglucinol/18-06_draft/Images/stained/ccoaomt_stained.png"))
fah1stained <- rasterGrob(readPNG("~/Documents/Uni/Phloroglucinol/18-06_draft/Images/stained/fah1_stained.png"))
omt1stained <- rasterGrob(readPNG("~/Documents/Uni/Phloroglucinol/18-06_draft/Images/stained/omt1_stained.png"))
ccr1stained <- rasterGrob(readPNG("~/Documents/Uni/Phloroglucinol/18-06_draft/Images/stained/ccr1_stained.png"))
cad4stained <- rasterGrob(readPNG("~/Documents/Uni/Phloroglucinol/18-06_draft/Images/stained/cad4_stained.png"))
cad5stained <- rasterGrob(readPNG("~/Documents/Uni/Phloroglucinol/18-06_draft/Images/stained/cad5_stained.png"))
cad4x5stained <- rasterGrob(readPNG("~/Documents/Uni/Phloroglucinol/18-06_draft/Images/stained/cad4x5_stained.png"))

WTunstained <- rasterGrob(readPNG("~/Documents/Uni/Phloroglucinol/18-06_draft/Images/unstained/WT_unstained.png"))
cl1unstained <- rasterGrob(readPNG("~/Documents/Uni/Phloroglucinol/18-06_draft/Images/unstained/4cl1_unstained.png"))
cl2unstained <- rasterGrob(readPNG("~/Documents/Uni/Phloroglucinol/18-06_draft/Images/unstained/4cl2_unstained.png"))
cl1x2unstained <- rasterGrob(readPNG("~/Documents/Uni/Phloroglucinol/18-06_draft/Images/unstained/4cl1x2_unstained.png"))
ccoaomtunstained <- rasterGrob(readPNG("~/Documents/Uni/Phloroglucinol/18-06_draft/Images/unstained/ccoaomt_unstained.png"))
fah1unstained <- rasterGrob(readPNG("~/Documents/Uni/Phloroglucinol/18-06_draft/Images/unstained/fah1_unstained.png"))
omt1unstained <- rasterGrob(readPNG("~/Documents/Uni/Phloroglucinol/18-06_draft/Images/unstained/omt1_unstained.png"))
ccr1unstained <- rasterGrob(readPNG("~/Documents/Uni/Phloroglucinol/18-06_draft/Images/unstained/ccr1_unstained.png"))
cad4unstained <- rasterGrob(readPNG("~/Documents/Uni/Phloroglucinol/18-06_draft/Images/unstained/cad4_unstained.png"))
cad5unstained <- rasterGrob(readPNG("~/Documents/Uni/Phloroglucinol/18-06_draft/Images/unstained/cad5_unstained.png"))
cad4x5unstained <- rasterGrob(readPNG("~/Documents/Uni/Phloroglucinol/18-06_draft/Images/unstained/cad4x5_unstained.png"))

WTfire <- rasterGrob(readPNG("~/Documents/Uni/Phloroglucinol/18-06_draft/Images/fire/WT_fire.png"))
cl1fire <- rasterGrob(readPNG("~/Documents/Uni/Phloroglucinol/18-06_draft/Images/fire/4cl1_fire.png"))
cl2fire <- rasterGrob(readPNG("~/Documents/Uni/Phloroglucinol/18-06_draft/Images/fire/4cl2_fire.png"))
cl1x2fire <- rasterGrob(readPNG("~/Documents/Uni/Phloroglucinol/18-06_draft/Images/fire/4cl1x2_fire.png"))
ccoaomtfire <- rasterGrob(readPNG("~/Documents/Uni/Phloroglucinol/18-06_draft/Images/fire/ccoaomt_fire.png"))
fah1fire <- rasterGrob(readPNG("~/Documents/Uni/Phloroglucinol/18-06_draft/Images/fire/fah1_fire.png"))
omt1fire <- rasterGrob(readPNG("~/Documents/Uni/Phloroglucinol/18-06_draft/Images/fire/omt1_fire.png"))
ccr1fire <- rasterGrob(readPNG("~/Documents/Uni/Phloroglucinol/18-06_draft/Images/fire/ccr1_fire.png"))
cad4fire <- rasterGrob(readPNG("~/Documents/Uni/Phloroglucinol/18-06_draft/Images/fire/cad4_fire.png"))
cad5fire <- rasterGrob(readPNG("~/Documents/Uni/Phloroglucinol/18-06_draft/Images/fire/cad5_fire.png"))
cad4x5fire <- rasterGrob(readPNG("~/Documents/Uni/Phloroglucinol/18-06_draft/Images/fire/cad4x5_fire.png"))
ccr1xfah1fire <- rasterGrob(readPNG("~/Documents/Uni/Phloroglucinol/18-06_draft/Images/fire/ccr1xfah1_fire.png"))
scale <- rasterGrob(readPNG("~/Documents/Uni/Phloroglucinol/18-06_draft/Images/fire/scale_aldehyde.png"))

# STAINED
pdf("stained_montage.pdf", height = 3.065, width = 10)
stained_pre <- plot_grid(
  cl1stained,
  cl2stained,
  cl1x2stained,
  ccoaomtstained,
  fah1stained,
  omt1stained,
  ccr1stained,
  cad4stained,
  cad5stained,
  cad4x5stained,
  labels = c("4cl1", "4cl2", "4cl1x4cl2", "ccoaomt1", "fah1", "omt1", "ccr1", "cad4", "cad5", "cad4xcad5"),
  label_fontfamily = "Helvetica",
  label_fontface = 3,
  scale = 0.98,
  ncol = 5,
  nrow = 2,
  hjust = 0,
  vjust = 1,
  label_x = 0.02,
  label_y = 0.98,
  label_size = 12)
stained <- plot_grid(
  WTstained,
  stained_pre,
  labels = c("Col-0", ""),
  label_fontfamily = "Helvetica",
  label_fontface = 1,
  scale = 0.99,
  ncol = 2,
  nrow = 1,
  hjust = 0.5,
  vjust = 1,
  label_x = 0.5,
  label_y = 0.98,
  label_size = 12,
  rel_widths = c(2.86, 7.14)
)
stained
dev.off()
system("gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/default -dNOPAUSE -dQUIET -dBATCH -dDetectDuplicateImages -dCompressFonts=true -r150 -sOutputFile=stained_mont.pdf stained_montage.pdf")


# UNSTAINED
pdf("unstained_montage.pdf", height = 3.065, width = 10)
unstained_pre <- plot_grid(
  cl1unstained,
  cl2unstained,
  cl1x2unstained,
  ccoaomtunstained,
  fah1unstained,
  omt1unstained,
  ccr1unstained,
  cad4unstained,
  cad5unstained,
  cad4x5unstained,
  labels = c("4cl1", "4cl2", "4cl1x4cl2", "ccoaomt1", "fah1", "omt1", "ccr1", "cad4", "cad5", "cad4xcad5"),
  label_fontfamily = "Helvetica",
  label_fontface = 3,
  scale = 0.98,
  ncol = 5,
  nrow = 2,
  hjust = 0,
  vjust = 1,
  label_x = 0.02,
  label_y = 0.98,
  label_size = 12)

unstained <- plot_grid(
  WTunstained,
  unstained_pre,
  labels = c("Col-0", ""),
  label_fontfamily = "Helvetica",
  label_fontface = 1,
  scale = 0.99,
  ncol = 2,
  nrow = 1,
  hjust = 0.5,
  vjust = 1,
  label_x = 0.5,
  label_y = 0.99,
  label_size = 12,
  rel_widths = c(2.86, 7.14)
)
unstained
dev.off()
system("gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/default -dNOPAUSE -dQUIET -dBATCH -dDetectDuplicateImages -dCompressFonts=true -r150 -sOutputFile=unstained_mont.pdf unstained_montage.pdf")

pdf("montage.pdf", height = 2, width = 10.3)
genotypes_mont <- plot_grid(
  stained,
  unstained,
  labels = c("",""),
  nrow = 2,
  ncol = 1
)
genotypes_mont
dev.off()
system("gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/default -dNOPAUSE -dQUIET -dBATCH -dDetectDuplicateImages -dCompressFonts=true -r150 -sOutputFile=mont.pdf montage.pdf")


# FIRE
pdf("fire_montage.pdf", height = 6, width = 8)
fire <- plot_grid(
  WTfire,
  cl1fire,
  cl2fire,
  cl1x2fire,
  ccoaomtfire,
  fah1fire,
  omt1fire,
  ccr1fire,
  ccr1xfah1fire,
  cad4fire,
  cad5fire,
  cad4x5fire,
  labels = c("Col-0", "4cl1", "4cl2", "4cl1x4cl2", "ccoaomt1", "fah1", "omt1", "ccr1", "ccr1xfah1", "cad4", "cad5", "cad4xcad5"),
  label_fontfamily = "Helvetica",
  label_fontface = 3,
  label_colour = "white",
  scale = 0.98,
  ncol = 4,
  nrow = 3,
  hjust = 0.5,
  vjust = 1,
  label_x = 0.5,
  label_y = 0.98)
fire
dev.off()
system("gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/printer -dNOPAUSE -dQUIET -dBATCH -dDetectDuplicateImages -dCompressFonts=true -r150 -sOutputFile=fire_mont.pdf fire_montage.pdf")



