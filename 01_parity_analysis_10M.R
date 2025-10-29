# ============================================================
# Batch analyzer for 10,000,000 digits of: π, e, φ, √2, √3, ln 2, ln 10, Catalan G, Euler γ
# ============================================================

mem.maxVSize(320000000000)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ---------------- Helpers ----------------
read_first_N_digits <- function(path, N = 1e7) {
  txt <- paste(readLines(path, warn = FALSE), collapse = "")
  digs <- gsub("\\D", "", txt)
  if (nchar(digs) < N) stop(sprintf("File %s has only %d digits", path, nchar(digs)))
  digsN <- substr(digs, 1, N)
  as.integer(strsplit(digsN, "", fixed = TRUE)[[1]])
}

parity_vec   <- function(d) as.integer(d %% 2)                  # 1=odd, 0=even
cum_dev      <- function(p) cumsum(ifelse(p == 1L, 1L, -1L))    # D(n)
runs_count   <- function(p) if (length(p)) 1L + sum(diff(p)!=0) else 0L

wilson_ci <- function(x, n, conf = 0.95) {
  if (n == 0) return(c(NA_real_, NA_real_))
  z <- qnorm(1 - (1 - conf)/2)
  ph <- x/n
  den <- 1 + z^2/n
  ctr <- (ph + z^2/(2*n))/den
  half <- z * sqrt((ph*(1-ph) + z^2/(4*n))/n)/den
  c(ctr - half, ctr + half)
}

transitions <- function(p) {
  if (length(p) < 2) return(setNames(rep(0L, 4), c("EE","EO","OE","OO")))
  a <- p[-length(p)]*2 + p[-1]
  c(EE = sum(a==0), EO = sum(a==1), OE = sum(a==2), OO = sum(a==3))
}

sliding_podd <- function(p, w = 1e4L) { # default 10k window for 10M digits
  n <- length(p); if (n < w) return(data.table(i = integer(), p_odd = numeric()))
  s <- c(0L, cumsum(p))
  k <- (w:n)
  podd <- (s[k+1] - s[k-w+1]) / w
  data.table(i = k, p_odd = podd)
}

block_props <- function(p, block = 1e6L) { # 10 blocks for 10M
  n <- length(p); B <- n %/% block
  if (B == 0) return(data.table(block = integer(), p_odd = numeric()))
  idx <- rep(1:B, each = block)
  idx <- idx[seq_len(B*block)]
  data.table(block = idx, odd = p[seq_len(B*block)])[
    , .(p_odd = mean(odd==1L)), by = block]
}

block_parity_test <- function(p, k = 2L) {
  n <- length(p); if (n < k) return(list(k=k, chisq=NA, df=NA, p=NA))
  mat <- embed(p, k)[, k:1, drop = FALSE]
  idx <- as.integer(mat %*% (2^(0:(k-1))))
  tab <- tabulate(idx + 1L, nbins = 2^k)
  E <- sum(tab)/(2^k)
  chisq <- sum((tab - E)^2 / E)
  list(k = k, chisq = chisq, df = (2^k - 1),
       p = pchisq(chisq, df = 2^k - 1, lower.tail = FALSE),
       counts = tab, expected = E)
}

digit_freq <- function(d) {
  tab <- tabulate(d + 1L, nbins = 10L)
  N <- length(d)
  exp <- rep(N/10, 10L)
  chisq <- sum((tab - exp)^2 / exp)
  pval <- pchisq(chisq, df = 9, lower.tail = FALSE)
  list(counts = tab, chisq = chisq, df = 9, p = pval)
}

monobit_test <- function(p) {
  n <- length(p); if (n == 0) return(list(z = NA_real_, p = NA_real_, S = NA_integer_))
  S <- sum(2L*p - 1L)
  z <- abs(S) / sqrt(n)
  pval <- 2 * (1 - pnorm(z))
  list(z = z, p = pval, S = S)
}

runs_test_nist <- function(p) {
  n <- length(p); if (n < 2) return(list(z = NA_real_, p = NA_real_, R = NA_integer_))
  phat <- mean(p)
  if (abs(phat - 0.5) > 0.5 - 2/sqrt(n)) return(list(z = NA_real_, p = NA_real_, R = NA_integer_))
  R <- 1L + sum(diff(p) != 0L)
  ER <- 1 + 2*n*phat*(1 - phat)
  VarR <- 2*n*phat*(1 - phat) * (2*n*phat*(1 - phat) - 1) / (n - 1)
  z <- (R - ER) / sqrt(VarR)
  pval <- 2 * (1 - pnorm(abs(z)))
  list(z = z, p = pval, R = R, ER = ER, VarR = VarR)
}

serial_test <- function(p, k = 2L) {
  n <- length(p); if (n < k) return(list(k = k, chisq = NA_real_, df = NA_integer_, p = NA_real_))
  mm <- embed(c(p, p[seq_len(k-1L)]), k)[, k:1, drop = FALSE]
  idx <- as.integer(mm %*% (2^(0:(k-1)))) + 1L
  tab <- tabulate(idx, nbins = 2^k)
  E <- length(idx) / (2^k)
  chisq <- sum((tab - E)^2 / E)
  df <- 2^k - 1
  list(k = k, chisq = chisq, df = df, p = pchisq(chisq, df = df, lower.tail = FALSE),
       counts = tab, expected = E)
}

approx_entropy <- function(p, m = 2L) {
  count_patterns <- function(pp, m) {
    n <- length(pp)
    mm <- embed(c(pp, pp[seq_len(m-1L)]), m)[, m:1, drop = FALSE]
    idx <- as.integer(mm %*% (2^(0:(m-1)))) + 1L
    tab <- tabulate(idx, nbins = 2^m)
    tab / n
  }
  Cm  <- count_patterns(p, m)
  Cmp <- count_patterns(p, m + 1L)
  phi_m  <- mean(log(pmax(Cm,  .Machine$double.eps)))
  phi_mp <- mean(log(pmax(Cmp, .Machine$double.eps)))
  list(m = m, ApEn = phi_m - phi_mp)
}

analyze_constant <- function(digits, name, win = 1e4L, block = 1e6L, kmax = 5L, sample_every = 1e5L) {
  stopifnot(all(digits %in% 0:9))
  p <- parity_vec(digits); N <- length(p)
  # Parity summary
  x_odd <- sum(p == 1L); p_odd <- x_odd / N; ci <- wilson_ci(x_odd, N)
  D <- cum_dev(p); scaled_max <- max(abs(D)) / sqrt(N)
  tr <- transitions(p); tr_tot <- sum(tr); tr_exp <- rep(tr_tot/4,4)
  tr_chi <- if (tr_tot>0) sum((tr - tr_exp)^2 / tr_exp) else NA_real_
  tr_p   <- if (is.na(tr_chi)) NA_real_ else pchisq(tr_chi, df=3, lower.tail=FALSE)
  rcount <- runs_count(p)
  # Sliding + blocks
  swin <- sliding_podd(p, w = win)
  btab <- block_props(p, block = block)
  # k-block tests
  blks <- lapply(2:kmax, function(k) block_parity_test(p, k = k))
  # Digit frequencies
  d0 <- digit_freq(digits)
  # Downsample D(n) for plotting
  keep <- seq(sample_every, N, by = sample_every)
  D_samp <- data.table(n = keep, D = D[keep], name = name)
  
  mono <- monobit_test(p)
  runz <- runs_test_nist(p)
  ser2 <- serial_test(p, k = 2L)
  ser3 <- serial_test(p, k = 3L)
  apen <- approx_entropy(p, m = 2L)
  
  summary <- data.table(
    name, N,
    p_odd = p_odd, ci_low = ci[1], ci_high = ci[2],
    runs = rcount,
    trans_EE = tr["EE"], trans_EO = tr["EO"], trans_OE = tr["OE"], trans_OO = tr["OO"],
    trans_chisq = tr_chi, trans_df = 3, trans_p = tr_p,
    max_abs_D = max(abs(D)), max_abs_D_over_sqrtN = scaled_max,
    digit_chisq = d0$chisq, digit_df = d0$df, digit_p = d0$p,
    mono_z = mono$z, mono_p = mono$p,
    runs_z = runz$z, runs_p = runz$p,
    serial2_chisq = ser2$chisq, serial2_df = ser2$df, serial2_p = ser2$p,
    serial3_chisq = ser3$chisq, serial3_df = ser3$df, serial3_p = ser3$p,
    apen_m2 = apen$ApEn
  )
  
  list(
    name = name, summary = summary,
    D_samples = D_samp, sliding = swin, blocks = btab,
    block_tests = blks, digit_counts = d0$counts
  )
}

# ---------------- Paths ----------------
paths <- list(
  "π"                 = "pi.txt",
  "e"                 = "e.txt",
  "φ"                 = "phi.txt",
  "√2"                = "sqrt2.txt",
  "√3"                = "sqrt3.txt",
  "ln 2"              = "ln2.txt",
  "ln 10"             = "ln10.txt",
  "Catalan G"         = "catalanG.txt",
  "Euler γ"           = "euler_gamma.txt",
  "Lemniscate \u03D6" = "lemniscate.txt"
)

# ---------------- Run all (10,000,000 digits each) ----------------
WIN <- 1e4L; BLOCK <- 1e6L; KMAX <- 5L; SAMPLE_EVERY <- 1e5L
results <- list()

for (nm in names(paths)) {
  message("Analyzing: ", nm)
  d <- read_first_N_digits(paths[[nm]], N = 1e7)   # first 10,000,000 digits
  results[[nm]] <- analyze_constant(d, name = nm, win = WIN, block = BLOCK, kmax = KMAX, sample_every = SAMPLE_EVERY)
}

# ---------------- Collect tables ----------------
all_summary <- rbindlist(lapply(results, `[[`, "summary"))
bp_table <- function(res) data.table(
  name = res$name,
  k    = sapply(res$block_tests, `[[`, "k"),
  chisq= sapply(res$block_tests, `[[`, "chisq"),
  df   = sapply(res$block_tests, `[[`, "df"),
  p    = sapply(res$block_tests, `[[`, "p")
)
all_k <- rbindlist(lapply(results, bp_table))
dig_table <- function(res) {
  dt <- data.table(digit = 0:9, count = as.integer(res$digit_counts))
  dt[, name := res$name]
  dt
}
all_digits <- rbindlist(lapply(results, dig_table))
fwrite(all_summary, "summary_10M.csv")
fwrite(all_k,       "kblock_tests_10M.csv")
fwrite(all_digits,  "digit_counts_10M.csv")

# ---------------- Plots (faceted) ----------------
# 1) Cumulative deviation D(n) sampled every 1e5 digits
D_all <- rbindlist(lapply(results, `[[`, "D_samples"))
p_cumdev <- ggplot(D_all, aes(n, D)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ name, scales = "free_y", ncol = 3) +
  labs(x = "n (digits)", y = "D(n) = #odd − #even",
       title = "Cumulative deviation D(n) across constants (10 million digits each)") +
  theme_classic()
ggsave("fig01.png", p_cumdev, width = 12, height = 9, dpi = 300)

# 2) Sliding-window odd share (w = 10,000)
summarize_slide <- function(dt, nbins = 5000L) {
  rng <- range(dt$i)
  brk <- seq(rng[1], rng[2], length.out = nbins + 1L)
  dt[, bin := cut(i, brk, include.lowest = TRUE)]
  dt[, .(
    i_mid = mean(i),
    m  = mean(p_odd),
    lo = quantile(p_odd, 0.05),
    hi = quantile(p_odd, 0.95)
  ), by = .(name, bin)]
}

S_all <- rbindlist(lapply(results, function(res) { tmp <- copy(res$sliding); tmp$name <- res$name; tmp }))
S_sum <- summarize_slide(S_all, nbins = 5000L)

p_slide_band <- ggplot(S_sum, aes(i_mid, m)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.25, fill = "gray70") +
  geom_line(linewidth = 0.1) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  facet_wrap(~ name, ncol = 3) +
  coord_cartesian(ylim = c(0.485, 0.515)) +
  labs(
    x = "n (digit index)",
    y = "Proportion odd (w = 10,000)",
    title = "Sliding-window odd share: mean ± 5–95% band"
  ) +
  theme_classic()
ggsave("fig02.png", p_slide_band, width = 12, height = 9, dpi = 300)

# 3) Digit-frequency bars 0–9 (proportion), faceted
dev <- copy(all_digits)[
  , prop := count / sum(count), by = name
][
  , `:=`(
    digit    = factor(digit, levels = 0:9),
    dev_prop = prop - 0.1,           
    dev_bp   = (prop - 0.1) * 1e4    
  )
]

# Auto-zoom around 0.1 so tiny bars are visible
max_bp <- max(abs(dev$dev_bp))
ylim <- 0.1 + c(-1, 1) * 1.2 * max_bp / 1e4  

p <- ggplot(dev, aes(digit)) +
  geom_segment(aes(xend = digit, y = 0.1, yend = 0.1 + dev_prop),
               linewidth = 6, colour = "red", lineend = "butt") +
  geom_hline(yintercept = 0.1, linetype = "solid", linewidth = .75) +
  facet_wrap(~ name, ncol = 3) +
  coord_cartesian(ylim = ylim) +
  scale_y_continuous(
    name = "Digit proportion",
    sec.axis = sec_axis(~ (. - 0.1) * 1e4, name = "Deviation (basis points)")
  ) +
  labs(x = "Digit",
       title = "Digit frequencies vs. uniform: red bars = deviation; dashed line = 0.1") +
  theme_classic() +
  theme(
    axis.title.y.right = element_text(color = "red"),
    axis.text.y.right  = element_text(color = "red"),
    axis.ticks.y.right = element_line(color = "red")
  )
ggsave("fig03.png", p, width = 12, height = 9, dpi = 300)

# ---------------- Abstract-ready blurbs ----------------
fmt <- function(x) ifelse(is.na(x), "NA", sprintf("%.6f", x))
cat("\nAbstract blurbs (10M each):\n")
for (r in split(all_summary, seq_len(nrow(all_summary)))) {
  r <- as.list(r)
  cat(sprintf(
    "%s: p_odd=%s (95%% CI %s–%s); max|D|/√N=%s; transitions χ²(3)=%.2f (p=%s); digits χ²(9)=%.2f (p=%s).\n",
    r$name, fmt(r$p_odd), fmt(r$ci_low), fmt(r$ci_high),
    fmt(r$max_abs_D_over_sqrtN), r$trans_chisq, fmt(r$trans_p),
    r$digit_chisq, fmt(r$digit_p)
  ))
}

cat("\nRandomness tests (10M each):\n")
for (r in split(all_summary, seq_len(nrow(all_summary)))) {
  r <- as.list(r)
  cat(sprintf(
    "%s: monobit z=%0.3f (p=%s); runs z=%0.3f (p=%s); serial2 χ²=%0.2f (p=%s); serial3 χ²=%0.2f (p=%s); ApEn(m=2)=%0.4f\n",
    r$name, r$mono_z, fmt(r$mono_p), r$runs_z, fmt(r$runs_p),
    r$serial2_chisq, fmt(r$serial2_p), r$serial3_chisq, fmt(r$serial3_p), r$apen_m2
  ))
}

