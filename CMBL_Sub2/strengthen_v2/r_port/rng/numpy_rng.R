## numpy_rng.R - numpy-compatible random generator in pure R (R_PORT_PLAN.md, gate G0)
##
## Reproduces numpy.random.default_rng (PCG64 seeded through SeedSequence) bit for bit for the subset of the Generator API that the
## version-2 pipeline uses: random, integers, standard_normal / normal, shuffle / permutation, choice (with and without replacement).
## Everything is plain R (doubles that hold exact integers; 128-bit and 64-bit values are kept as 16-bit and 32-bit limbs).
## Usage:   g <- np_rng(2026)            # default_rng(2026);  np_rng(c(2026, 5, 1, 7)) is default_rng(SeedSequence([2026, 5, 1, 7]))
##          np_random(g, 5); np_integers(g, 0, 17, 5); np_standard_normal(g, 9); np_permutation(g, 12); np_choice_noreplace_idx(g, 20, 6)
## The generator object is an environment: it is modified in place, like numpy's.
## Tests: test_numpy_rng.R (against vectors written by make_numpy_reference.py from numpy 2.4.4).
## Speed note: %% and %/% on doubles are slow in R, so the hot paths use floor(x / 2^k) instead (exact for the magnitudes used here).

.P32 <- 4294967296
.lib_dir <- (function() {       # folder of this file also when it is sourced from another sourced file (the innermost source() frame whose ofile is this file)
  for (fr in rev(sys.frames())) { f <- fr$ofile; if (is.character(f) && grepl("numpy_rng[.]R$", f)) return(dirname(normalizePath(f))) }
  f <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(f)) dirname(normalizePath(f)) else getwd()
})()

## ---- 32-bit helpers (unsigned values held as doubles) ---------------------------------------------------------------
.u32_mul <- function(a, b) {
  ah <- a %/% 65536; al <- a %% 65536; bh <- b %/% 65536; bl <- b %% 65536
  (((ah * bl + al * bh) %% 65536) * 65536 + al * bl) %% .P32
}
.u32_xor <- function(a, b) {
  bitwXor(as.integer(a %/% 65536), as.integer(b %/% 65536)) * 65536 + bitwXor(as.integer(a %% 65536), as.integer(b %% 65536))
}
.u32_sub <- function(a, b) (a - b) %% .P32

## ---- SeedSequence (numpy/random/bit_generator.pyx) ---------------------------------------------------------------------
.seed_words <- function(seed) {
  out <- numeric(0)
  for (v in as.numeric(seed)) {
    if (v < 0 || v != floor(v)) stop("seed values must be non-negative integers")
    if (v == 0) out <- c(out, 0) else while (v > 0) { out <- c(out, v %% .P32); v <- v %/% .P32 }
  }
  out
}
.hashmix <- function(value, hc) {
  value <- .u32_xor(value, hc)
  hc <- .u32_mul(hc, 0x931e8875)
  value <- .u32_mul(value, hc)
  value <- .u32_xor(value, value %/% 65536)
  c(value, hc)
}
.mix <- function(x, y) {
  r <- .u32_sub(.u32_mul(0xca01f9dd, x), .u32_mul(0x4973f715, y))
  .u32_xor(r, r %/% 65536)
}
np_seedseq_pool <- function(entropy) {
  ent <- .seed_words(entropy); pool <- numeric(4); hc <- 0x43b0d7e5
  for (i in 1:4) { hm <- .hashmix(if (i <= length(ent)) ent[i] else 0, hc); pool[i] <- hm[1]; hc <- hm[2] }
  for (src in 1:4) for (dst in 1:4) if (src != dst) {
    hm <- .hashmix(pool[src], hc); hc <- hm[2]; pool[dst] <- .mix(pool[dst], hm[1])
  }
  if (length(ent) > 4) for (src in 5:length(ent)) for (dst in 1:4) {
    hm <- .hashmix(ent[src], hc); hc <- hm[2]; pool[dst] <- .mix(pool[dst], hm[1])
  }
  pool
}
## generate_state(n64, uint64) as 32-bit words: word 2k-1 is the low and word 2k the high half of the k-th 64-bit value
np_seedseq_words <- function(pool, n64 = 4) {
  n <- 2 * n64; w <- numeric(n); hc <- 0x8b51f9dd
  for (i in 1:n) {
    dv <- pool[((i - 1) %% 4) + 1]; dv <- .u32_xor(dv, hc); hc <- .u32_mul(hc, 0x58f38ded)
    dv <- .u32_mul(dv, hc); dv <- .u32_xor(dv, dv %/% 65536); w[i] <- dv
  }
  w
}

## ---- 128-bit scalars as 8 little-endian 16-bit limbs -----------------------------------------------------------------
.M16 <- c(0xF645, 0x9FCC, 0xDF64, 0x4385, 0x5DA4, 0x1FC6, 0xED05, 0x2360)   # PCG_DEFAULT_MULTIPLIER_128 = 0x2360ED051FC65DA44385DF649FCCF645
.carry128 <- function(r) {
  cr <- 0
  for (k in 1:8) { t <- r[k] + cr; cr <- floor(t / 65536); r[k] <- t - cr * 65536 }
  r
}
.add128 <- function(a, b) .carry128(a + b)
.mul128 <- function(a, b) {
  r <- numeric(8)
  for (i in 1:8) if (a[i] != 0) r[i:8] <- r[i:8] + a[i] * b[1:(9 - i)]
  .carry128(r)
}
.step_scalar <- function(s, inc) .add128(.mul128(s, .M16), inc)
.shl1_or1 <- function(x) { cr <- 0; for (k in 1:8) { t <- x[k] * 2 + cr; cr <- floor(t / 65536); x[k] <- t - cr * 65536 }; x[1] <- x[1] + 1; x }

## ---- vectorised PCG64 (XSL-RR 128/64): many lanes in lock step, jump-ahead by affine maps ------------------------------
.toep <- function(a) { m <- matrix(0, 8, 8); for (i in 1:8) for (k in i:8) m[i, k] <- a[k - i + 1]; m }
.carry_mat <- function(R) {
  for (k in 1:7) { v <- R[, k]; cr <- floor(v / 65536); R[, k] <- v - cr * 65536; R[, k + 1] <- R[, k + 1] + cr }
  v <- R[, 8]; R[, 8] <- v - floor(v / 65536) * 65536
  R
}
.affine_compose <- function(f, g) list(a = .mul128(f$a, g$a), c = .add128(.mul128(f$a, g$c), f$c))   # f after g
.affine_pow <- function(f, n) {
  r <- list(a = c(1, 0, 0, 0, 0, 0, 0, 0), c = numeric(8)); b <- f
  while (n > 0) { if (n %% 2 == 1) r <- .affine_compose(b, r); n <- n %/% 2; if (n > 0) b <- .affine_compose(b, b) }
  r
}
.apply_affine <- function(S, f) .carry_mat(S %*% .toep(f$a) + matrix(f$c, nrow(S), 8, byrow = TRUE))
.POW2 <- 2^(0:15)
## output function on states (n x 8 limbs): rotr64(hi64 xor lo64, state >> 122); returns the 64-bit value as hi32 and lo32
.pcg_output <- function(S) {
  n <- nrow(S); rows <- seq_len(n)
  Xv <- c(bitwXor(as.integer(S[, 1]), as.integer(S[, 5])) + 0, bitwXor(as.integer(S[, 2]), as.integer(S[, 6])) + 0,
          bitwXor(as.integer(S[, 3]), as.integer(S[, 7])) + 0, bitwXor(as.integer(S[, 4]), as.integer(S[, 8])) + 0)
  rot <- floor(S[, 8] / 1024); q <- floor(rot / 16); p <- .POW2[rot - q * 16 + 1]; ip <- 65536 / p
  yq <- vector("list", 4); yr <- vector("list", 4)
  for (j in 1:4) {
    idx <- j + q; idx <- idx - 4 * (idx > 4)
    y <- Xv[rows + (idx - 1) * n]; a <- floor(y / p)
    yq[[j]] <- a; yr[[j]] <- y - a * p
  }
  o1 <- yq[[1]] + yr[[2]] * ip; o2 <- yq[[2]] + yr[[3]] * ip; o3 <- yq[[3]] + yr[[4]] * ip; o4 <- yq[[4]] + yr[[1]] * ip
  list(lo = o1 + o2 * 65536, hi = o3 + o4 * 65536)
}
## B lanes (a power of two) x L steps: returns B*L outputs in stream order and the state after them
.pcg_block <- function(state, inc, B, L) {
  f1 <- list(a = .M16, c = inc)
  S <- matrix(state, 1, 8); lanes <- 1; fj <- .affine_pow(f1, L)
  while (lanes < B) { S <- rbind(S, .apply_affine(S, fj)); lanes <- lanes * 2; fj <- .affine_compose(fj, fj) }
  Mt <- .toep(.M16); Cm <- matrix(inc, B, 8, byrow = TRUE)
  HI <- matrix(0, B, L); LO <- matrix(0, B, L)
  for (l in 1:L) { S <- .carry_mat(S %*% Mt + Cm); o <- .pcg_output(S); HI[, l] <- o$hi; LO[, l] <- o$lo }
  list(hi = as.vector(t(HI)), lo = as.vector(t(LO)), state = S[B, ])
}

## ---- generator object ----------------------------------------------------------------------------------------------------
np_rng <- function(seed) {
  g <- new.env(parent = emptyenv())
  w <- np_seedseq_words(np_seedseq_pool(seed))
  limbs <- function(lo32, hi32) c(lo32 %% 65536, lo32 %/% 65536, hi32 %% 65536, hi32 %/% 65536)
  s0 <- c(limbs(w[3], w[4]), limbs(w[1], w[2]))       # state = val[0] * 2^64 + val[1]
  i0 <- c(limbs(w[7], w[8]), limbs(w[5], w[6]))       # inc   = val[2] * 2^64 + val[3]
  inc <- .shl1_or1(i0)
  st <- .step_scalar(numeric(8), inc)                 # pcg_setseq_128_srandom_r: state = 0; step; state += initstate; step
  st <- .step_scalar(.add128(st, s0), inc)
  g$state <- st; g$inc <- inc
  g$hi <- numeric(0); g$lo <- numeric(0); g$pos <- 1; g$n <- 0; g$blk <- 4096; g$epoch <- 0
  g$has32 <- FALSE; g$u32 <- 0
  g
}
np_state <- function(g) list(state = g$state, inc = g$inc)
.ensure <- function(g, n) {
  avail <- g$n - g$pos + 1
  while (avail < n) {
    K <- max(g$blk, n - avail)
    B <- min(2048, 2^ceiling(log2(sqrt(K)))); L <- ceiling(K / B)
    blk <- .pcg_block(g$state, g$inc, B, L)
    g$state <- blk$state
    keep <- g$n >= g$pos
    g$hi <- c(if (keep) g$hi[g$pos:g$n] else numeric(0), blk$hi)
    g$lo <- c(if (keep) g$lo[g$pos:g$n] else numeric(0), blk$lo)
    g$pos <- 1; g$n <- length(g$hi); avail <- g$n; g$epoch <- g$epoch + 1
    g$blk <- min(1048576, g$blk * 4)
  }
  invisible()
}
np_raw64 <- function(g, n) { .ensure(g, n); p <- g$pos; idx <- p:(p + n - 1); g$pos <- p + n; list(hi = g$hi[idx], lo = g$lo[idx]) }

## ---- 32-bit stream (numpy's pcg64_next32 keeps the unused high half for the next 32-bit draw) -------------------------------
.next32 <- function(g) {
  if (g$has32) { g$has32 <- FALSE; return(g$u32) }
  if (g$pos > g$n) .ensure(g, 1)
  p <- g$pos; g$pos <- p + 1; g$has32 <- TRUE; g$u32 <- g$hi[p]
  g$lo[p]
}
.peek32 <- function(g, m) {
  out <- numeric(m); i <- 0
  if (g$has32) { out[1] <- g$u32; i <- 1 }
  rem <- m - i
  if (rem > 0) {
    nw <- (rem + 1) %/% 2; .ensure(g, nw); p <- g$pos; idx <- p:(p + nw - 1)
    v <- numeric(2 * nw); v[c(TRUE, FALSE)] <- g$lo[idx]; v[c(FALSE, TRUE)] <- g$hi[idx]
    out[(i + 1):m] <- v[seq_len(rem)]
  }
  out
}
.consume32 <- function(g, d) {
  if (d <= 0) return(invisible())
  s <- if (g$has32) 1 else 0
  if (d <= s) { g$has32 <- FALSE; return(invisible()) }
  rem <- d - s; nw <- (rem + 1) %/% 2; .ensure(g, nw)
  last <- g$pos + nw - 1; g$pos <- g$pos + nw
  if (rem %% 2 == 1) { g$has32 <- TRUE; g$u32 <- g$hi[last] } else g$has32 <- FALSE
  invisible()
}
.next_double <- function(g) {
  if (g$pos > g$n) .ensure(g, 1)
  p <- g$pos; g$pos <- p + 1
  (g$hi[p] * 2097152 + floor(g$lo[p] / 2048)) * 2^-53
}

## ---- random(), integers() ------------------------------------------------------------------------------------------------
np_random <- function(g, n) { w <- np_raw64(g, n); (w$hi * 2097152 + floor(w$lo / 2048)) * 2^-53 }
.mul_split <- function(u, r) {              # u * r (u < 2^32, r <= 2^32) as high and low 32 bits
  uh <- floor(u / 65536); ul <- u - uh * 65536; rh <- floor(r / 65536); rl <- r - rh * 65536
  t1 <- ul * rl; t2 <- uh * rl + ul * rh; t3 <- uh * rh
  t2h <- floor(t2 / 65536); low <- t1 + (t2 - t2h * 65536) * 65536; lh <- floor(low / .P32)
  list(hi = t3 + t2h + lh, lo = low - lh * .P32)
}
.lemire_scalar <- function(g, rng) {        # one bounded draw in [0, rng], rng >= 1 (numpy buffered_bounded_lemire_uint32)
  rexcl <- rng + 1; thr <- (4294967295 - rng) %% rexcl
  repeat { m <- .mul_split(.next32(g), rexcl); if (m$lo >= thr) return(m$hi) }
}
np_integers <- function(g, low, high, n) {  # Generator.integers(low, high, n), int64, endpoint = FALSE
  rng <- high - low - 1
  if (rng < 0) stop("high must be above low")
  if (rng == 0) return(rep(low, n))
  if (rng > 4294967295) stop("ranges above 2^32 - 1 are not implemented")
  if (rng == 4294967295) { out <- numeric(n); for (i in seq_len(n)) out[i] <- .next32(g); return(low + out) }
  rexcl <- rng + 1; thr <- (4294967295 - rng) %% rexcl
  res <- numeric(n); got <- 0
  while (got < n) {
    need <- n - got
    m <- .mul_split(.peek32(g, need), rexcl)
    acc <- m$lo >= thr
    .consume32(g, need)
    k <- sum(acc)
    if (k > 0) res[(got + 1):(got + k)] <- m$hi[acc]
    got <- got + k
  }
  low + res
}

## ---- shuffle, permutation, choice ----------------------------------------------------------------------------------------
.random_interval <- function(g, mx) {       # numpy random_interval: masked rejection, value in [0, mx], mx >= 1
  mask1 <- 2^ceiling(log2(mx + 1))
  repeat { v <- .next32(g) %% mask1; if (v <= mx) return(v) }
}
np_shuffle <- function(g, x, first = 1) {   # Generator.shuffle on a 1-D vector (returns the shuffled copy); first = 1 shuffles everything
  n <- length(x)
  if (n > first) for (i in (n - 1):first) { j <- .random_interval(g, i); if (j != i) { t <- x[i + 1]; x[i + 1] <- x[j + 1]; x[j + 1] <- t } }
  x
}
np_permutation <- function(g, n) np_shuffle(g, 0:(n - 1))          # Generator.permutation(n), 0-based values
np_permute <- function(g, x) np_shuffle(g, x)                      # Generator.permutation(array)
.shuffle_int <- function(g, x, first = 1) {  # numpy _shuffle_int (used by choice): bounded Lemire draws, not random_interval
  n <- length(x)
  if (n > first) for (i in (n - 1):first) { j <- .lemire_scalar(g, i); t <- x[i + 1]; x[i + 1] <- x[j + 1]; x[j + 1] <- t }
  x
}
np_choice_noreplace_idx <- function(g, pop, k, shuffle = TRUE) {   # choice(arange(pop), k, replace = False): 0-based indices
  if (k > pop) stop("cannot take a larger sample than population when replace is False")
  if (k == 0) return(numeric(0))
  if (pop > 10000 && k > pop %/% (if (shuffle) 50 else 20)) {       # tail shuffle
    idx <- .shuffle_int(g, 0:(pop - 1), first = max(pop - k, 1))
    return(idx[(pop - k + 1):pop])
  }
  idx <- numeric(k); used <- logical(pop); base <- pop - k          # Floyd's algorithm
  for (j in base:(pop - 1)) {
    v <- if (j == 0) 0 else .lemire_scalar(g, j)
    if (!used[v + 1]) { used[v + 1] <- TRUE; idx[j - base + 1] <- v } else { used[j + 1] <- TRUE; idx[j - base + 1] <- j }
  }
  if (shuffle && k > 1) idx <- .shuffle_int(g, idx, 1)
  idx
}
np_choice_replace_idx <- function(g, pop, n) np_integers(g, 0, pop, n)

## ---- standard_normal (ziggurat, 52-bit; tables from numpy's ziggurat_constants.h) -------------------------------------------
.zig <- new.env(parent = emptyenv())
np_load_ziggurat <- function(path = file.path(.lib_dir, "ziggurat_nor.f64")) {
  if (!file.exists(path)) stop("ziggurat table file not found: ", path)
  v <- readBin(path, "double", n = 770, size = 8)
  .zig$ki <- v[1:256]; .zig$wi <- v[257:512]; .zig$fi <- v[513:768]; .zig$r <- v[769]; .zig$inv_r <- v[770]
  invisible(TRUE)
}
.zig_get <- function() { if (is.null(.zig$ki)) np_load_ziggurat(); .zig }
.zdecode <- function(hi, lo) {
  l8 <- floor(lo / 256); l9 <- floor(lo / 512); h29 <- floor(hi / 536870912)
  list(idx = lo - l8 * 256, sgn = l8 - floor(l8 / 2) * 2, rabs = (hi - h29 * 536870912) * 8388608 + l9)
}
.norm_slow <- function(g, idx, sgn, rabs, x) {   # the rest of numpy's loop for a draw that failed the fast test
  z <- .zig_get()
  repeat {
    if (idx == 0) {
      repeat {
        xx <- -z$inv_r * log1p(-.next_double(g)); yy <- -log1p(-.next_double(g))
        if (yy + yy > xx * xx) return(if ((rabs %/% 256) %% 2 == 1) -(z$r + xx) else z$r + xx)
      }
    } else if ((z$fi[idx] - z$fi[idx + 1]) * .next_double(g) + z$fi[idx + 1] < exp(-0.5 * x * x)) return(x)
    w <- np_raw64(g, 1); d <- .zdecode(w$hi, w$lo)
    idx <- d$idx; sgn <- d$sgn; rabs <- d$rabs
    x <- rabs * z$wi[idx + 1]; if (sgn == 1) x <- -x
    if (rabs < z$ki[idx + 1]) return(x)
  }
}
np_standard_normal <- function(g, n) {
  z <- .zig_get(); out <- numeric(n); got <- 0
  while (got < n) {
    need <- n - got
    m <- need + ceiling(need * 0.03) + 64                   # decode a little more than needed: slow draws consume extra words
    .ensure(g, m + 256)
    p <- g$pos; ep <- g$epoch; iw <- p:(p + m - 1)
    d <- .zdecode(g$hi[iw], g$lo[iw])
    x <- d$rabs * z$wi[d$idx + 1] * (1 - 2 * d$sgn)
    slowv <- which(!(d$rabs < z$ki[d$idx + 1])); ns <- length(slowv)
    q <- 1; si <- 1; rebased <- FALSE
    while (got < n && q <= m) {
      while (si <= ns && slowv[si] < q) si <- si + 1
      nxt <- if (si <= ns) slowv[si] else m + 1
      cnt <- min(nxt - q, n - got)
      if (cnt > 0) { out[(got + 1):(got + cnt)] <- x[q:(q + cnt - 1)]; got <- got + cnt; q <- q + cnt }
      if (got >= n || q > m) break
      g$pos <- p + q                                         # word q is a slow draw: consume it, then run the rest of the loop
      out[got + 1] <- .norm_slow(g, d$idx[q], d$sgn[q], d$rabs[q], x[q]); got <- got + 1
      if (g$epoch != ep) { rebased <- TRUE; break }          # the buffer was refilled: chunk coordinates are void
      q <- g$pos - p + 1
    }
    if (!rebased) g$pos <- p + q - 1
  }
  out
}
np_normal <- function(g, loc, scale, n) loc + scale * np_standard_normal(g, n)
