## pyrandom.R - Python's `random.Random(seed)` (Mersenne Twister MT19937 as in CPython: init_by_array seeding for an int, getrandbits, _randbelow_with_getrandbits, shuffle) in pure R,
## for the selection scripts that shuffle with random.Random(2026).shuffle (select_datasets.py, select_heldout.py, select_rnaseq.py). 32-bit words are doubles; products are done in 16-bit limbs.
## py_random(seed) -> generator (an environment); py_getrandbits(g, k) k <= 32; py_randbelow(g, n); py_shuffle(g, x) returns the shuffled copy (as random.shuffle does in place).
.M32 <- 4294967296
.pr_hi <- function(x) floor(x / 65536); .pr_lo <- function(x) x - floor(x / 65536) * 65536
.pr_xor <- function(a, b) bitwXor(.pr_hi(a), .pr_hi(b)) * 65536 + bitwXor(.pr_lo(a), .pr_lo(b))
.pr_and <- function(a, b) bitwAnd(.pr_hi(a), .pr_hi(b)) * 65536 + bitwAnd(.pr_lo(a), .pr_lo(b))
.pr_shr <- function(x, s) floor(x / 2^s)
.pr_shl <- function(x, s) { y <- x * 2^s; y - floor(y / .M32) * .M32 }
.pr_mul <- function(a, b) { ah <- .pr_hi(a); al <- .pr_lo(a); bh <- .pr_hi(b); bl <- .pr_lo(b); cross <- (ah * bl + al * bh) %% 65536; (al * bl + cross * 65536) %% .M32 }        # a * b mod 2^32
.pr_init_genrand <- function(s) { mt <- numeric(624); mt[1] <- s; for (i in 2:624) mt[i] <- (.pr_mul(1812433253, .pr_xor(mt[i - 1], .pr_shr(mt[i - 1], 30))) + (i - 1)) %% .M32; mt }
.pr_init_by_array <- function(key) {
  mt <- .pr_init_genrand(19650218); N <- 624L; i <- 2L; j <- 1L; klen <- length(key)             # i, j are 1-based positions of C's i (1) and j (0)
  for (k in seq_len(max(N, klen))) {
    mt[i] <- (.pr_xor(mt[i], .pr_mul(.pr_xor(mt[i - 1], .pr_shr(mt[i - 1], 30)), 1664525)) + key[j] + (j - 1)) %% .M32
    i <- i + 1L; j <- j + 1L; if (i > N) { mt[1] <- mt[N]; i <- 2L }; if (j > klen) j <- 1L
  }
  for (k in seq_len(N - 1L)) {
    mt[i] <- (.pr_xor(mt[i], .pr_mul(.pr_xor(mt[i - 1], .pr_shr(mt[i - 1], 30)), 1566083941)) - (i - 1)) %% .M32
    i <- i + 1L; if (i > N) { mt[1] <- mt[N]; i <- 2L }
  }
  mt[1] <- 2147483648; mt
}
py_random <- function(seed) {                                                                       # random.Random(int seed): the key is the 32-bit little-endian words of abs(seed)
  n <- abs(seed); key <- numeric(0); repeat { key <- c(key, n %% .M32); n <- floor(n / .M32); if (n == 0) break }
  g <- new.env(); g$mt <- .pr_init_by_array(key); g$mti <- 625L; g                                   # mti > 624: generate on first use
}
.pr_twist <- function(g) {
  mt <- g$mt; N <- 624L; M <- 397L; mag <- c(0, 2567483615)
  for (kk in 1:(N - M)) { y <- (if (mt[kk] >= 2147483648) 2147483648 else 0) + (mt[kk + 1] %% 2147483648); mt[kk] <- .pr_xor(.pr_xor(mt[kk + M], .pr_shr(y, 1)), mag[(y %% 2) + 1]) }
  for (kk in (N - M + 1):(N - 1)) { y <- (if (mt[kk] >= 2147483648) 2147483648 else 0) + (mt[kk + 1] %% 2147483648); mt[kk] <- .pr_xor(.pr_xor(mt[kk + M - N], .pr_shr(y, 1)), mag[(y %% 2) + 1]) }
  y <- (if (mt[N] >= 2147483648) 2147483648 else 0) + (mt[1] %% 2147483648); mt[N] <- .pr_xor(.pr_xor(mt[M], .pr_shr(y, 1)), mag[(y %% 2) + 1])
  g$mt <- mt; g$mti <- 1L
}
.pr_next32 <- function(g) {
  if (g$mti > 624L) .pr_twist(g)
  y <- g$mt[g$mti]; g$mti <- g$mti + 1L
  y <- .pr_xor(y, .pr_shr(y, 11)); y <- .pr_xor(y, .pr_and(.pr_shl(y, 7), 2636928640)); y <- .pr_xor(y, .pr_and(.pr_shl(y, 15), 4022730752)); .pr_xor(y, .pr_shr(y, 18))
}
py_getrandbits <- function(g, k) { stopifnot(k >= 1, k <= 32); .pr_shr(.pr_next32(g), 32 - k) }
py_randbelow <- function(g, n) { if (n < 1) return(0); k <- if (n == 1) 1 else floor(log2(n)) + 1; r <- py_getrandbits(g, k); while (r >= n) r <- py_getrandbits(g, k); r }       # n.bit_length() = floor(log2(n)) + 1 for n >= 1
py_shuffle <- function(g, x) { for (i in rev(seq_len(length(x) - 1L))) { j <- py_randbelow(g, i + 1) + 1; tmp <- x[i + 1]; x[i + 1] <- x[j]; x[j] <- tmp }; x }
