"""patch_r_copies.py - prepares the two scratch copies of Fig_5abcd_Gene_Interaction_Network_Analysis.r for running with Rscript:
  * font registration (the script asks for family = "Arial"; the Windows graphics device of Rscript only knows it after windowsFonts(); RStudio does this implicitly),
  * svg saves skipped (the composite figure uses the png panels only).
The patch is inserted after the line that prints 'Packages loaded'. It is idempotent. The scientific code of the scripts is not touched (the alias copy differs from the as-reported copy only by the Slc2a2 -> Glut2 mapping of the STRING node names)."""
import os
HERE = os.path.dirname(os.path.abspath(__file__))
PATCH = '''# ---- WP10 copy: environment patch (8 Nov 2026) ----
# Rscript opens the pdf device for plots that the script prints (a bare plot object on its own line); the pdf device does not know family = "Arial" (invalid font type), RStudio uses a device that does: any implicitly opened device is ragg here.
options(device = function(...) ragg::agg_png(filename = tempfile(fileext = ".png"), width = 8, height = 7.5, units = "in", res = 72))
if (.Platform$OS.type == "windows") windowsFonts(Arial = windowsFont("Arial"))
.orig_ggsave <- ggplot2::ggsave
ggsave <- function(filename, ...) { if (grepl("[.]svg$", filename)) { cat("  (svg save skipped in this copy)\\n"); return(invisible(NULL)) }; .orig_ggsave(filename, ...) }
# ---- end of the environment patch ----
'''
for v in ('as_reported', 'alias'):
    p = os.path.join(HERE, v, 'Fig_5abcd_Gene_Interaction_Network_Analysis.r')
    lines = open(p, encoding='utf-8').read().split('\n')
    if any('environment patch' in l for l in lines):          # replace an older version of the patch block
        a = [i for i, l in enumerate(lines) if '---- WP10 copy: environment patch' in l][0]
        b = [i for i, l in enumerate(lines) if '---- end of the environment patch ----' in l][0]
        del lines[a:b + 1]
    idx = [i for i, l in enumerate(lines) if 'Packages loaded' in l]
    assert len(idx) == 1, (v, idx)
    lines[idx[0] + 1:idx[0] + 1] = PATCH.rstrip('\n').split('\n')
    open(p, 'w', encoding='utf-8', newline='\n').write('\n'.join(lines))
    print(v, 'patched after line', idx[0] + 1)
