#!/usr/bin/env Rscript

required_packages <- c("rEEMSplots", "maps")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required R package is not installed: ", pkg)
  }
}

library(rEEMSplots)
library(maps)

project_dir <- "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/chamaea_wgs"
eems_dir <- file.path(project_dir, "4_analyses", "05_eems")
mcmcpath <- file.path(eems_dir, "chamaea_eems_run1")
plot_dir <- file.path(eems_dir, "plots")

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(mcmcpath)) {
  stop("EEMS MCMC directory does not exist: ", mcmcpath)
}

eems.plots(
  mcmcpath = mcmcpath,
  plotpath = file.path(plot_dir, "chamaea_eems_standard"),
  longlat = TRUE,
  out.png = FALSE,
  add.outline = TRUE,
  add.grid = TRUE,
  add.demes = TRUE,
  res = 600
)

eems.plots(
  mcmcpath = mcmcpath,
  plotpath = file.path(plot_dir, "chamaea_eems_state-outlines"),
  longlat = TRUE,
  out.png = FALSE,
  add.outline = TRUE,
  add.grid = TRUE,
  add.demes = TRUE,
  res = 600,
  col.demes = "black",
  col.grid = "gray20",
  m.plot.xy = {
    map(
      "state",
      regions = c("california", "oregon", "nevada"),
      col = "black",
      lwd = 2,
      add = TRUE
    )
  },
  q.plot.xy = {
    map(
      "state",
      regions = c("california", "oregon", "nevada"),
      col = "black",
      lwd = 2,
      add = TRUE
    )
  }
)

eems.plots(
  mcmcpath = mcmcpath,
  plotpath = file.path(plot_dir, "chamaea_eems_state-outlines-notitle"),
  longlat = TRUE,
  out.png = FALSE,
  add.outline = TRUE,
  add.grid = TRUE,
  add.demes = TRUE,
  add.title = FALSE,
  res = 600,
  col.demes = "black",
  col.grid = "gray20",
  m.plot.xy = {
    map(
      "state",
      regions = c("california", "oregon", "nevada"),
      col = "black",
      lwd = 2,
      add = TRUE
    )
  },
  q.plot.xy = {
    map(
      "state",
      regions = c("california", "oregon", "nevada"),
      col = "black",
      lwd = 2,
      add = TRUE
    )
  }
)
