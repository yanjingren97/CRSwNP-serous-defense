options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({library(grid); library(ragg)})

out <- Sys.getenv("ARTICLE_WORKFLOW_CLEAN_OUT", unset = file.path(".", "deliverables", "Figure1_article_style_clean"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "QA_final_size_96dpi"), recursive = TRUE, showWarnings = FALSE)

font <- "Arial"
ink <- "#25323D"
muted <- "#687986"
line <- "#B8C5CD"
navy <- "#315B78"
teal <- "#4D9C91"
blue <- "#4C8FB3"
orange <- "#D88755"
gold <- "#C9A44D"
panel <- "#F7F9FA"
blue_light <- "#EAF3F8"
teal_light <- "#EAF5F3"
orange_light <- "#FBF0E9"
gold_light <- "#F8F3E6"
grey_light <- "#F1F4F5"

rr <- function(x, y, w, h, fill = "white", col = line, lwd = .8,
               lty = "solid", r = .012) {
  grid.roundrect(unit(x, "npc"), unit(y, "npc"), unit(w, "npc"), unit(h, "npc"),
                 r = unit(r, "npc"),
                 gp = gpar(fill = fill, col = col, lwd = lwd, lty = lty))
}

txt <- function(label, x, y, size = 8, col = ink, face = "plain",
                just = "centre", lineheight = 1.06) {
  grid.text(label, unit(x, "npc"), unit(y, "npc"), just = just,
            gp = gpar(fontfamily = font, fontsize = size, col = col,
                      fontface = face, lineheight = lineheight))
}

header <- function(label, x, y, w, colour, size = 8.3) {
  grid.roundrect(unit(x, "npc"), unit(y, "npc"), unit(w, "npc"), unit(.095, "npc"),
                 r = unit(.011, "npc"), gp = gpar(fill = colour, col = colour))
  txt(label, x, y, size, "white", "bold")
}

subbox <- function(label, x, y, w, h, fill, col, size = 7.4, dashed = FALSE,
                   face = "plain") {
  rr(x, y, w, h, fill, col, .7, if (dashed) "22" else "solid", .010)
  txt(label, x, y, size, ink, face, lineheight = 1.04)
}

arrow_between <- function(x1, x2, y = .605, dashed = FALSE) {
  grid.lines(unit(c(x1, x2), "npc"), unit(c(y, y), "npc"),
             arrow = arrow(type = "closed", length = unit(1.8, "mm")),
             gp = gpar(col = muted, lwd = 1.0,
                       lty = if (dashed) "22" else "solid", lineend = "round"))
}

draw <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))

  xs <- c(.085, .245, .430, .675, .885)
  ws <- c(.140, .140, .200, .240, .170)
  card_y <- .615
  card_h <- .650
  for (i in seq_along(xs)) rr(xs[i], card_y, ws[i], card_h, panel, line, .85, "solid", .012)

  header("PUBLIC COHORTS", xs[1], .895, ws[1], navy, 7.7)
  subbox("GSE235711\npaired design\npolyp–ethmoid", xs[1], .720, .116, .135,
         blue_light, blue, 7.5)
  subbox("GSE276503\nexternal\nmulti-site cohort", xs[1], .535, .116, .135,
         blue_light, blue, 7.5)
  subbox("BIOLOGICAL\nUNIT\ndonor or biopsy", xs[1], .365, .116, .095,
         grey_light, muted, 6.7, FALSE, "bold")

  header("FIXED PROCESSING", xs[2], .895, ws[2], teal, 7.1)
  txt("Fixed QC\nand annotation", xs[2], .755, 7.7, ink, "bold")
  grid.lines(unit(c(xs[2] - .052, xs[2] + .052), "npc"), unit(c(.710, .710), "npc"),
             gp = gpar(col = line, lwd = .7))
  txt("Disease-\nindependent\nsecretory states", xs[2], .640, 7.4, ink)
  txt("Integration\nreciprocal label\ntransfer", xs[2], .520, 7.5, ink)
  txt("Pre-specified\nexclusions\nunit-level inference", xs[2], .395, 6.8, muted)

  header("DISCOVERY", xs[3], .895, ws[3], blue, 8.5)
  subbox("COMPOSITION\nANALYSIS\nSerous-glandular\nfractions\nState remodeling\nUnit-level contrasts",
         xs[3], .700, .170, .175, blue_light, blue, 7.15)
  rr(xs[3], .495, .170, .190, orange_light, orange, .7, "solid", .010)
  txt("WITHIN-CELL\nTRANSCRIPTIONAL\nANALYSIS",
      xs[3], .540, 6.9, ink, "plain", lineheight = 1.00)
  txt("Serous pseudobulk\nDefinition genes excluded\nCross-cohort concordance",
      xs[3], .450, 6.3, ink, "plain", lineheight = 1.00)
  subbox("LOCKED 40-GENE\nSEROUS-DEFENSE\nMODULE",
         xs[3], .335, .170, .085, orange_light, orange, 6.9, FALSE, "bold")

  header("CONTEXT & VALIDATION", xs[4], .895, ws[4], orange, 7.8)
  txt("FROM THE LOCKED MODULE", xs[4], .815, 6.9, muted, "bold")
  subbox("FUNCTIONAL &\nREGULATORY CONTEXT\nGO / Reactome\nState coordinate\nCollecTRI / PROGENy",
         xs[4], .705, .205, .145, teal_light, teal, 7.05)
  subbox("EXTERNAL BULK\nVALIDATION\nGSE36830 · GSE136825\nLocked score\nPaired analysis",
         xs[4], .520, .205, .145, gold_light, gold, 7.05)
  subbox("CROSS-COHORT\nCO-EXPRESSION\nGSE136825 → GSE36830\nReproducibility",
         xs[4], .370, .205, .125, blue_light, blue, 6.9)

  header("BOUNDARY AUDITS", xs[5], .895, ws[5], gold, 8.0)
  txt("EXPLORATORY BRANCH", xs[5], .815, 6.9, muted, "bold")
  subbox("Bulk-cohort\nECM coupling", xs[5], .700, .142, .105,
         grey_light, muted, 7.3, TRUE)
  subbox("GSE235714\nGeoMx proxy", xs[5], .535, .142, .105,
         grey_light, muted, 7.3, TRUE)
  subbox("PXD013330\nsecretion\nproteomics", xs[5], .370, .142, .105,
         grey_light, muted, 7.3, TRUE)
  txt("Exploratory only", xs[5], .305, 6.9, muted, "bold")

  # Only the first three columns are a sequential processing chain.
  arrow_between(xs[1] + ws[1] / 2 + .004, xs[2] - ws[2] / 2 - .004)
  arrow_between(xs[2] + ws[2] / 2 + .004, xs[3] - ws[3] / 2 - .004)

  # Discovery, context/validation and boundary evidence converge rather than form a false sequence.
  for (x in xs[3:5]) {
    grid.lines(unit(c(x, x), "npc"), unit(c(.288, .250), "npc"),
               gp = gpar(col = if (x == xs[5]) gold else line, lwd = .85,
                         lty = if (x == xs[5]) "22" else "solid"))
  }
  grid.lines(unit(c(xs[3], xs[5]), "npc"), unit(c(.250, .250), "npc"),
             gp = gpar(col = line, lwd = .9))
  grid.lines(unit(c(.690, .690), "npc"), unit(c(.250, .215), "npc"),
             arrow = arrow(type = "closed", length = unit(1.7, "mm")),
             gp = gpar(col = muted, lwd = 1.0))

  txt("ANATOMY-AWARE EVIDENCE SYNTHESIS", .690, .190, 8.0, muted, "bold")
  subbox("VARIABLE REDUCTION IN\nSEROUS-GLANDULAR\nREPRESENTATION",
         .350, .105, .305, .095, blue_light, blue, 7.8, FALSE, "bold")
  subbox("REPRODUCIBLE\nWITHIN-CELL SUPPRESSION\nof the serous host-defense program",
         .680, .105, .305, .095, orange_light, orange, 7.8, FALSE, "bold")
  txt("Transcriptomic tissue-state program  |  no causal or treatment-response inference",
      .515, .035, 6.9, muted)
}

width_mm <- 180
height_mm <- 105

ragg::agg_tiff(file.path(out, "Figure1_article_style_clean_180mm.tiff"),
               width = width_mm, height = height_mm, units = "mm", res = 600,
               compression = "lzw", background = "white")
draw(); dev.off()

ragg::agg_png(file.path(out, "QA_final_size_96dpi", "Figure1_article_style_clean_QA96.png"),
              width = width_mm, height = height_mm, units = "mm", res = 96,
              background = "white")
draw(); dev.off()

cairo_pdf(file.path(out, "Figure1_article_style_clean_180mm.pdf"),
          width = width_mm / 25.4, height = height_mm / 25.4,
          family = font, bg = "white")
draw(); dev.off()

svg(file.path(out, "Figure1_article_style_clean_180mm.svg"),
    width = width_mm / 25.4, height = height_mm / 25.4,
    family = font, bg = "white")
draw(); dev.off()

file.copy(file.path("code", "plot_article_workflow_clean_final.R"),
          file.path(out, "plot_article_workflow_clean_final.R"), overwrite = TRUE)

writeLines(c(
  "Figure 1 clean five-column article workflow",
  "180 x 105 mm; Arial; 600 dpi LZW RGB TIFF.",
  "The visual style follows the supplied article_workflow.tiff reference.",
  "Only cohorts -> fixed processing -> discovery are shown as a sequential chain.",
  "Discovery, context/validation and boundary audits converge at the evidence-synthesis line.",
  "",
  capture.output(sessionInfo())
), file.path(out, "R_sessionInfo_and_design_notes.txt"))

message(normalizePath(out, winslash = "/", mustWork = TRUE))
