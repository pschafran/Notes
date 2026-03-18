# Patch for GENESPACE v1.3.1 run_mcscanx() bug:
# When MCScanX produces 0 syntenic blocks, the collinearity file contains no
# anchor lines. The grep in fread() then exits code 1, throwing an unhandled
# error. This patch wraps the fread() call in tryCatch() so that zero-block
# pairs return an empty result instead of crashing.
#
# Usage: source("patch_run_mcscanx.R") before calling run_genespace()

run_mcscanx_patched <- function(hits, blkSize, nGaps, tmpDir, MCScanX_hCall) {
    mcsID1 <- mcsID2 <- NULL
    hitCols <- c("ofID1", "ofID2", "chr1", "chr2", "ord1", "ord2", "score")
    if (!all(hitCols %in% colnames(hits)))
        stop("Looks like the hits data.table is misformed. It must be a data.table",
            "comparing two genomes with the columns:", hitCols, "\n")
    if (length(blkSize) != 1)
        stop("blkSize must be a integer of length 1\n")
    if (length(nGaps) != 1)
        stop("blkSize must be a integer of length 1\n")
    if (!is.integer(blkSize)) {
        blkSize <- as.integer(blkSize)
        if (is.na(blkSize))
            stop("Cannot coerce blkSize to an integer\n")
    }
    if (!is.integer(nGaps)) {
        nGaps <- as.integer(nGaps)
        if (is.na(nGaps))
            stop("Cannot coerce nGaps to an integer\n")
    }
    if (!file.exists(MCScanX_hCall))
        stop("Cannot find MCScanX_h executable in", MCScanX_hCall, "\n")
    tmpd <- file.path(tmpDir, paste0("tmp_", paste(sample(c(letters,
        LETTERS), 20, replace = T), collapse = "")))
    if (dir.exists(tmpd))
        unlink(tmpd, recursive = T)
    dir.create(tmpd)
    on.exit(expr = unlink(tmpd, recursive = T))
    ord1 <- ord2 <- chr1 <- chr2 <- ofID1 <- ofID2 <- NULL
    u1 <- subset(hits, !duplicated(ofID1))
    u2 <- subset(hits, !duplicated(ofID2))
    setkey(u1, chr1, ord1)
    setkey(u2, chr2, ord2)
    u1[, `:=`(mcsID1, paste0("aa", as.numeric(factor(chr1, levels = unique(chr1)))))]
    u2[, `:=`(mcsID2, paste0("bb", as.numeric(factor(chr2, levels = unique(chr2)))))]
    mcs1 <- u1[, c("mcsID1", "ofID1", "ord1", "ord1")]
    setnames(mcs1, c("chr", "id", "start", "end"))
    mcs2 <- u2[, c("mcsID2", "ofID2", "ord2", "ord2")]
    setnames(mcs2, c("chr", "id", "start", "end"))
    mcsGffIn <- rbind(mcs1, mcs2)
    mcsBlsIn <- hits[, c("ofID1", "ofID2", "score")]
    mcsBlsIn[, `:=`(ofID2, paste0(ofID2, "xxxx"))]
    mcsGffIn$id[grepl("bb", mcsGffIn$chr)] <- paste0(mcsGffIn$id[grepl("bb",
        mcsGffIn$chr)], "xxxx")
    blFile <- file.path(tmpd, "mcs.homology")
    gfFile <- file.path(tmpd, "mcs.gff")
    colFile <- file.path(tmpd, "mcs.collinearity")
    fwrite(mcsGffIn, file = gfFile, sep = "\t", quote = FALSE,
        col.names = FALSE, showProgress = FALSE, verbose = FALSE)
    fwrite(mcsBlsIn, file = blFile, sep = "\t", quote = FALSE,
        col.names = FALSE, showProgress = FALSE, verbose = FALSE)
    mcsCom <- sprintf("-a -b 2 -c 2 -m %s -s %s %s", nGaps, blkSize,
        file.path(tmpd, "mcs"))
    comout <- system2(MCScanX_hCall, mcsCom, stdout = TRUE, stderr = TRUE)
    idg <- strsplit(as.character(hits$ofID1[1]), "_")[[1]][1]

    # PATCHED: wrap fread in tryCatch to handle zero-block case
    # When MCScanX finds no syntenic blocks, grep exits code 1 -> fread errors
    collin <- tryCatch(
        suppressWarnings(fread(
            cmd = sprintf("cat %s | grep %s_ | grep :", colFile, idg),
            col.names = c("blkID", "gn1", "gn2"),
            select = 1:3, showProgress = FALSE, header = FALSE
        )),
        error = function(e) {
            data.table(blkID = character(0), gn1 = character(0), gn2 = character(0))
        }
    )

    if (nrow(collin) > 1) {
        blkID <- gn1 <- gn2 <- NULL
        collin[, `:=`(blkID, as.numeric(sapply(blkID, function(x) strsplit(x,
            "-")[[1]][1])) + 1)]
        collin[, `:=`(gn2, gsub("xxxx$", "", gn2))]
        mcsb <- collin$blkID
        names(mcsb) <- with(collin, paste(gn1, gn2))
        return(mcsb)
    }
}

# Apply the patch
environment(run_mcscanx_patched) <- asNamespace("GENESPACE")
assignInNamespace("run_mcscanx", run_mcscanx_patched, ns = "GENESPACE")

cat("GENESPACE::run_mcscanx patched: zero-block pairs will now be handled gracefully.\n")
