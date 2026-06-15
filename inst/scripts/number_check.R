# Extraction of all numerical values cited in the manuscript and the
# response-to-editor file, taken directly from the regenerated .rds data.
# Run from the rcpmcp package root. Output is printed to the console and
# written to inst/scripts/number_check_output.txt.

data_main <- "inst/scripts/data"
data_si   <- "inst/sup_info/data"
out_path  <- "inst/scripts/number_check_output.txt"
if (file.exists(out_path)) file.remove(out_path)

emit <- function(...) {
  line <- paste0(...)
  cat(line, "\n", sep = "")
  cat(line, "\n", sep = "", file = out_path, append = TRUE)
}
r3   <- function(x) formatC(round(as.numeric(x), 3), format = "f", digits = 3)
r4   <- function(x) formatC(round(as.numeric(x), 4), format = "f", digits = 4)
vec3 <- function(x) paste(r3(x), collapse = ", ")
rng3 <- function(x) paste0("[", r3(min(x)), ", ", r3(max(x)), "]")

gap_info <- function(adj, unadj) {
  g <- unadj - adj
  k <- which.max(g)
  paste0("adjusted=", rng3(adj), "  unadjusted=", rng3(unadj),
         "  max gap=", r3(max(g)), " at K=", k)
}

emit("==================================================================")
emit(" Manuscript / Response number extraction (regenerated .rds)")
emit("==================================================================")

## ================= Figure 1 / Section 3.2 (rho = 0) =================
f1d <- readRDS(file.path(data_main, "fig1_data.rds"))
emit("")
emit("################ Figure 1 / Section 3.2 (rho = 0) ################")
emit("N (K=1..5): ", paste(round(f1d$fig1_dat_01$N), collapse = ", "))
emit("(N is identical for both f1 values: ",
     identical(round(f1d$fig1_dat_01$N), round(f1d$fig1_dat_02$N)), ")")
for (tg in c("fig1_dat_01", "fig1_dat_02")) {
  d <- f1d[[tg]]
  emit("--- ", d$f1, " ---")
  emit("  disjunctive power, formula (K=1..5): ", vec3(d$Power_F))
  emit("  null RCP M1 formula (K=1..5): ", vec3(d$RCP0_M1_F),
       "  | K=1=", r3(d$RCP0_M1_F[1]), " -> K=5=", r3(d$RCP0_M1_F[5]))
  emit("  null RCP M2 formula (K=1..5): ", vec3(d$RCP0_M2_F),
       "  | K=1=", r3(d$RCP0_M2_F[1]), " -> K=5=", r3(d$RCP0_M2_F[5]))
}

## ================= Figure 2 / Section 3.3 (rho = 0) =================
f2 <- readRDS(file.path(data_main, "fig2_data.rds"))$fig2adj_df
get_panel <- function(panel, f1, approach = NULL) {
  s <- f2[f2$panel == panel & f2$f1 == f1, ]
  if (!is.null(approach)) s <- s[!is.na(s$Approach) & s$Approach == approach, ]
  s <- s[order(s$K), ]
  s$value
}
emit("")
emit("################ Figure 2 / Section 3.3 (rho = 0) ################")
for (f1 in c("f[1]==0.1", "f[1]==0.2")) {
  emit("--- ", f1, " ---")
  gm1 <- get_panel("gamma[M1]~(adjusted)", f1)
  gm2 <- get_panel("gamma[M2]~(adjusted)", f1)
  emit("  gamma_M1,K (K=1..5): ", vec3(gm1), "  | K=1=", r3(gm1[1]), " -> K=5=", r3(gm1[5]))
  emit("  gamma_M2,K (K=1..5): ", vec3(gm2), "  | K=1=", r3(gm2[1]), " -> K=5=", r3(gm2[5]))
}
# Reference values (single-endpoint null RCP, K = 1) from Figure 1 formula data
emit("Reference values phi_Ref (= null RCP at K=1, formula):")
emit("  f1=0.1: phi_Ref^(1)=", r3(f1d$fig1_dat_01$RCP0_M1_F[1]),
     "  phi_Ref^(2)=", r3(f1d$fig1_dat_01$RCP0_M2_F[1]))
emit("  f1=0.2: phi_Ref^(1)=", r3(f1d$fig1_dat_02$RCP0_M1_F[1]),
     "  phi_Ref^(2)=", r3(f1d$fig1_dat_02$RCP0_M2_F[1]))
# Bottom row: RCP under H1, adjusted vs unadjusted (from Figure 1 H1 data, rho=0)
emit("RCP under H1, adjusted vs unadjusted (bottom row of Figure 2):")
emit("  M1, f1=0.1: ", gap_info(f1d$fig1_dat_01$RCP1_M1_Adj, f1d$fig1_dat_01$RCP1_M1_Unadj))
emit("  M1, f1=0.2: ", gap_info(f1d$fig1_dat_02$RCP1_M1_Adj, f1d$fig1_dat_02$RCP1_M1_Unadj))
emit("  M2, f1=0.1: ", gap_info(f1d$fig1_dat_01$RCP1_M2_Adj, f1d$fig1_dat_01$RCP1_M2_Unadj))
emit("  M2, f1=0.2: ", gap_info(f1d$fig1_dat_02$RCP1_M2_Adj, f1d$fig1_dat_02$RCP1_M2_Unadj))

## ================= Table 1 / Section 4 =================
t1 <- readRDS(file.path(data_main, "table1_data.rds"))
emit("")
emit("################ Table 1 / Section 4 ################")
emit("Parameters: delta=", paste(t1$APP_DELTA, collapse = ","),
     "  rho=", t1$APP_RHO, "  fs=", paste(t1$APP_FS, collapse = ","),
     "  alpha=", t1$APP_ALPHA)
emit("N (K=1..", t1$APP_K, "): ", paste(round(t1$N), collapse = ", "))
emit("gamma_M1,K: ", vec3(t1$gamma_M1_adj))
emit("gamma_M2,K: ", vec3(t1$gamma_M2_adj))
emit("null RCP M1 unadjusted: ", vec3(t1$null_unadj_M1))
emit("null RCP M2 unadjusted: ", vec3(t1$null_unadj_M2))
emit("null RCP M1 adjusted:   ", vec3(t1$null_adj_M1))
emit("null RCP M2 adjusted:   ", vec3(t1$null_adj_M2))
emit("H1 RCP M1 unadjusted:   ", vec3(t1$unadj_M1))
emit("H1 RCP M2 unadjusted:   ", vec3(t1$unadj_M2))
emit("H1 RCP M1 adjusted:     ", vec3(t1$adj_M1))
emit("H1 RCP M2 adjusted:     ", vec3(t1$adj_M2))
emit("Max H1 reduction (unadj - adj): M1=", r3(max(t1$unadj_M1 - t1$adj_M1)),
     "  M2=", r3(max(t1$unadj_M2 - t1$adj_M2)))

## ================= Figure 3 (RCP vs rho) =================
f3 <- readRDS(file.path(data_main, "fig3_data.rds"))
d3 <- f3$fig3rho_df
emit("")
emit("################ Figure 3 (RCP vs rho) ################")
emit("rho sequence: ", paste(f3$RHO_SEQ, collapse = ", "))
emit("K values present: ", paste(sort(unique(d3$K)), collapse = ", "))
for (mth in c("Method~1", "Method~2")) {
  for (kk in sort(unique(d3$K))) {
    s <- d3[d3$panel == mth & d3$K == kk & d3$Approach == "Formula", ]
    s <- s[order(s$rho), ]
    emit("  ", mth, ", K=", kk, ", formula RCP vs rho(", paste(s$rho, collapse=","), "): ", vec3(s$value))
  }
}

## ================= Figure 4 (null RCP under 4 MCPs, adjusted) =================
f4 <- readRDS(file.path(data_main, "fig4_data.rds"))
d4 <- f4$fig4_df
emit("")
emit("################ Figure 4 (null RCP, 4 MCPs, adjusted thresholds) ################")
emit("Reference lines phi_Ref:")
emit("  f1=0.1: M1=", r3(f4$fig4_ref$M1_01), "  M2=", r3(f4$fig4_ref$M2_01))
emit("  f1=0.2: M1=", r3(f4$fig4_ref$M1_02), "  M2=", r3(f4$fig4_ref$M2_02))
emit("Range of null RCP across K and MCP (should stay near phi_Ref):")
for (f1 in c("f[1]==0.1", "f[1]==0.2")) {
  s <- d4[d4$f1 == f1, ]
  emit("  ", f1, ": M1 ", rng3(s$RCP0_M1), "  M2 ", rng3(s$RCP0_M2))
  emit("    max |RCP0_M1 - phi_Ref|=",
       r3(max(abs(s$RCP0_M1 - if (f1 == "f[1]==0.1") f4$fig4_ref$M1_01 else f4$fig4_ref$M1_02))),
       "  max |RCP0_M2 - phi_Ref|=",
       r3(max(abs(s$RCP0_M2 - if (f1 == "f[1]==0.1") f4$fig4_ref$M2_01 else f4$fig4_ref$M2_02))))
}

## ================= Supporting Information S1-S6 =================
emit("")
emit("################ Supporting Information ################")

## S1: unknown variance (known vs unknown)
s1 <- readRDS(file.path(data_si, "siS1_data.rds"))$siS1_df
emit("--- S1 (unknown variance): max |Known - Unknown| ---")
for (hh in c("H0", "H1")) {
  for (mth in c("Method~1", "Method~2")) {
    kn <- s1[s1$H == hh & s1$panel == mth & s1$variance == "Known", ]
    un <- s1[s1$H == hh & s1$panel == mth & s1$variance == "Unknown", ]
    kn <- kn[order(kn$K), ]; un <- un[order(un$K), ]
    emit("  ", hh, " ", mth, ": max diff=", r4(max(abs(kn$value - un$value))))
  }
}

## S2: region-of-interest heterogeneity (r1 = 0, 0.1, 0.2)
s2 <- readRDS(file.path(data_si, "siS2_data.rds"))$siS2_df
emit("--- S2 (region-1 effect r1): RCP at K=1 and K=5 ---")
for (mth in c("Method~1", "Method~2")) {
  for (r1v in sort(unique(s2$r1))) {
    s <- s2[s2$panel == mth & s2$r1 == r1v, ]
    s <- s[order(s$K), ]
    emit("  ", mth, ", r1=", r1v, ": K=1=", r3(s$value[1]), " K=5=", r3(s$value[length(s$value)]))
  }
}

## S3: mixed endpoints (Homogeneous vs Mixed), formula vs simulation
s3 <- readRDS(file.path(data_si, "siS3_data.rds"))$siS3_df
emit("--- S3 (mixed endpoints): RCP range and formula-vs-sim max diff ---")
for (mth in c("Method~1", "Method~2")) {
  for (cf in c("Homogeneous", "Mixed")) {
    ff <- s3[s3$panel == mth & s3$config == cf & s3$Approach == "Formula", ]
    sm <- s3[s3$panel == mth & s3$config == cf & s3$Approach == "Simulation", ]
    ff <- ff[order(ff$K), ]; sm <- sm[order(sm$K), ]
    emit("  ", mth, " ", cf, ": formula ", rng3(ff$value),
         "  max|F-S|=", r4(max(abs(ff$value - sm$value))))
  }
}

## S4: correlation x null region (rho in 0,0.3,0.5,0.8)
s4 <- readRDS(file.path(data_si, "siS4_data.rds"))$siS4_df
emit("--- S4 (correlation x null region): RCP at K=5 by rho ---")
for (mth in c("Method~1", "Method~2")) {
  for (rv in sort(unique(s4$rho))) {
    s <- s4[s4$panel == mth & s4$rho == rv, ]
    s <- s[order(s$K), ]
    emit("  ", mth, ", rho=", rv, ": K=1=", r3(s$value[1]), " K=5=", r3(s$value[length(s$value)]))
  }
}

## S5: unequal variances (Equal vs Unequal), formula vs simulation
s5 <- readRDS(file.path(data_si, "siS5_data.rds"))$siS5_df
emit("--- S5 (unequal variances): H0 calibration and formula-vs-sim max diff ---")
for (hh in c("H0", "H1")) {
  for (mth in c("Method~1", "Method~2")) {
    for (st in c("Equal", "Unequal")) {
      ff <- s5[s5$H == hh & s5$panel == mth & s5$struct == st & s5$Approach == "Formula", ]
      sm <- s5[s5$H == hh & s5$panel == mth & s5$struct == st & s5$Approach == "Simulation", ]
      ff <- ff[order(ff$K), ]; sm <- sm[order(sm$K), ]
      emit("  ", hh, " ", mth, " ", st, ": formula ", rng3(ff$value),
           "  max|F-S|=", r4(max(abs(ff$value - sm$value))))
    }
  }
}

## S6: non-standard correlation (Independence / CS / AR(1)), formula vs simulation
s6 <- readRDS(file.path(data_si, "siS6_data.rds"))$siS6_df
emit("--- S6 (non-standard correlation): H1 RCP range and formula-vs-sim max diff ---")
for (hh in c("H0", "H1")) {
  for (mth in c("Method~1", "Method~2")) {
    for (st in sort(unique(s6$struct))) {
      ff <- s6[s6$H == hh & s6$panel == mth & s6$struct == st & s6$Approach == "Formula", ]
      sm <- s6[s6$H == hh & s6$panel == mth & s6$struct == st & s6$Approach == "Simulation", ]
      ff <- ff[order(ff$K), ]; sm <- sm[order(sm$K), ]
      emit("  ", hh, " ", mth, " ", st, ": formula ", rng3(ff$value),
           "  max|F-S|=", r4(max(abs(ff$value - sm$value))))
    }
  }
}

emit("")
emit("==================================================================")
emit(" Session information")
emit("==================================================================")
si <- capture.output(sessionInfo())
for (l in si) emit(l)
