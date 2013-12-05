addClusters.semirnd <-
function (ks, clusterStack, ratios, no.overlap = T, clust.size = 30, 
    allowed.genes = NA, dont.allow = NA) 
{
    cat.new("Generating ", length(ks), "additional semi-random clusters\n")
    cat.new("Adding them at indices:", ks, "\n")
    if (is.null(clusterStack$k)) 
        old.k <- 0
    else old.k <- clusterStack$k
    clusterStack$k <- NULL
    names(clusterStack) <- NULL
    new.clusts <- list()
    seed.genes <- NULL
    for (i in ks) {
        if (no.overlap) {
            if (length(clusterStack) > 0) {
                clusterStack$k <- length(clusterStack)
                seed.genes <- genes.not.in.any.clusters(clusterStack, 
                  rownames(ratios))
                clusterStack$k <- NULL
                names(clusterStack) <- NULL
            }
        }
        if (!is.na(allowed.genes)) 
            seed.genes <- allowed.genes
        else seed.genes <- rownames(ratios)
        if (!is.na(dont.allow) && length(dont.allow) > 0 && length(seed.genes) > 
            0) 
            seed.genes <- seed.genes[!seed.genes %in% dont.allow]
        if (length(seed.genes) < clust.size * 10) {
            clusterStack$k <- length(clusterStack)
            cnt <- genes.in.how.many.clusters(clusterStack, rownames(ratios))
            clusterStack$k <- NULL
            names(clusterStack) <- NULL
            srt <- names(cnt)[sort(cnt, index = T)$ix]
            seed.genes <- unique(c(seed.genes, srt[1:(clust.size * 
                10)]))
        }
        seed.gene <- NULL
        if (no.overlap) {
            non.genes <- seed.genes
            if (length(non.genes) > 0) 
                seed.gene <- sample(non.genes, 1)
        }
        clust <- gen.clust.semirandom(ratios, clust.size = clust.size, 
            ratios.thresh = ratios.thresh, seed.genes = seed.genes, 
            seed.row = seed.gene)
        names(clusterStack) <- NULL
        clust$k <- i
        clust$resid <- residOneClust(ratios[clust$rows, clust$cols], 
            varNorm = TRUE)
        clust$p.clust <- NA
        clust$net.p.old <- get.cluster.subnet.p.values(clust, 
            networks, gene.ids, edge.centric = T, net.q0 = net.q0)
        cat.new("Cluster", clust$k, ": nrow=", clust$nrows, "; ncol=", 
            clust$ncols, "; resid=", clust$resid, "\n")
        clusterStack[[i]] <- clust
        names(clusterStack) <- NULL
        if (no.overlap && !is.na(dont.allow) && length(dont.allow) > 
            0) 
            dont.allow <- unique(c(dont.allow, clust$rows))
        else dont.allow <- clust$rows
    }
    clusterStack$k <- length(clusterStack)
    return(clusterStack)
}
all.iter.params <-
structure(list(iter = 0, resid.scaling = 1, bad.move.temperature = 0.15, 
    motif.scaling = 0, min.unique.seqs = 3, max.unique.seqs = 200, 
    n.motifs = 1, min.motif.width = 9, max.motif.width = 25, 
    seq.length = 250, bg.order = 3, e.val.cutoff = 999999, p.val.cutoff = 0.9, 
    motif.e.val.cutoff = 999999, meme.model = structure("zoops", .Names = "1"), 
    dynamic.operon.shift = FALSE, use.mcast = FALSE, max.motif.gap = 200, 
    net.scaling = 0.1, min.motif.width.pal = 10, max.motif.width.pal = 25, 
    r0 = 1, p0 = 0, q0 = 0.1, net.q0 = structure(c(0.2, 0.08, 
    0, 0.01, 0, 0, 0, 0, 0, 0, 0.01), .Names = c("operons", "met", 
    "prolinks.RS", "prolinks.PP", "prolinks.GN", "prolinks.GC", 
    "COG.code", "cond.sims", "predictome.CP", "predictome.GF", 
    "predictome.PP"))), .Names = c("iter", "resid.scaling", "bad.move.temperature", 
"motif.scaling", "min.unique.seqs", "max.unique.seqs", "n.motifs", 
"min.motif.width", "max.motif.width", "seq.length", "bg.order", 
"e.val.cutoff", "p.val.cutoff", "motif.e.val.cutoff", "meme.model", 
"dynamic.operon.shift", "use.mcast", "max.motif.gap", "net.scaling", 
"min.motif.width.pal", "max.motif.width.pal", "r0", "p0", "q0", 
"net.q0"))
biclust.version <-
"1.1"
cat.new <-
function (..., log = get.global("out.log"), file = "", sep = " ", 
    fill = FALSE, labels = NULL, append = FALSE) 
{
    cat(..., file = file, sep = sep, fill = fill, labels = labels, 
        append = append)
    if (file == "" && !is.null(log) && log != "") 
        try(cat(..., file = log, append = T))
}
choose.node <-
function (node.names, node.coords) 
{
    cat.new("Select a node:\n")
    node <- node.names[identify(node.coords, n = 1, plot = F)]
    return(node)
}
CircleGraphLayout <-
function (nodes, radius = 100) 
{
    node.count <- length(nodes)
    r <- radius
    phi <- 2 * pi/node.count
    node.coords <- matrix(0, nrow = node.count, ncol = 2)
    node.coords[, 1] <- r * sin(1:node.count * phi)
    node.coords[, 2] <- r * cos(1:node.count * phi)
    invisible(node.coords)
}
cluster.correlation <-
function (clust1, clust2, ratios = get.global("ratios"), use.cols = c("union", 
    "intersection")) 
{
    if (use.cols == "intersection") 
        cols <- clust1$cols[clust1$cols %in% clust2$cols]
    else if (use.cols == "union") 
        cols <- c(clust1$cols, clust2$cols)
    rat1 <- ratios[clust1$rows, cols]
    prof1 <- apply(rat1, 2, mean, na.rm = T)
    rat2 <- ratios[clust2$rows, cols]
    prof2 <- apply(rat2, 2, mean, na.rm = T)
    cor <- cor(prof1, prof2, use = "pairwise")
    return(cor)
}
compare.cluster.motifs <-
function (cluster1, cluster2, rev.comp = T, e.v.cutoff = 10) 
{
    m.out1 <- cluster1$motif.out
    m.out2 <- cluster2$motif.out
    if (is.null(m.out1) || is.null(m.out2) || is.null(m.out1$pssms) || 
        is.null(m.out2$pssms) || length(m.out1$pssms) <= 0 || 
        length(m.out2$pssms) <= 0) 
        return(-1)
    cors <- matrix(-1, nrow = length(m.out1$pssms), ncol = length(m.out2$pssms))
    for (i in 1:length(m.out1$pssms)) {
        if (i > ncol(cors)) 
            break
        for (j in i:length(m.out2$pssms)) {
            if (!is.na(m.out1$e.values) && !is.na(m.out2$e.values) && 
                !is.na(m.out1$e.values[i]) && !is.na(m.out2$e.values[j]) && 
                m.out1$e.values[i] <= e.v.cutoff && m.out2$e.values[j] <= 
                e.v.cutoff) {
                cors[i, j] <- compare.pssms(m.out1$pssms[[i]], 
                  m.out2$pssms[[j]], rev.comp = rev.comp)
            }
        }
    }
    cors[cors == 0] <- -1
    return(cors)
}
condition.sim.pvalue <-
function (conds, networks, all.conds, edge.centric = T) 
{
    p.value <- subnet.p.value(conds, networks$cond.sims, all.conds, 
        edge.centric, n.genes = networks$n.conds, adj.list = networks$adj.lists$cond.sims)
    return(p.value)
}
date.biclust.run <-
"Wed Jun 14 10:20:47 2006"
edge.list.to.adj.list <-
function (edges) 
{
    out <- list()
    genes <- unique(as.vector(edges))
    for (gene in genes) {
        tmp <- unique(c(edges[, 2][edges[, 1] == gene], edges[, 
            1][edges[, 2] == gene]))
        out[[gene]] <- tmp
    }
    return(out)
}
edge.list.to.netmatrix <-
function (edges, edge.names = NULL) 
{
    if (is.null(edge.names)) 
        edge.names <- unique(c(edges[, 1], edges[, 2]))
    out <- matrix(0, nrow = length(edge.names), ncol = length(edge.names))
    rownames(out) <- colnames(out) <- edge.names
    for (i in edge.names) {
        joined <- edges[, 2][which(edges[, 1] == i)]
        joined <- joined[which(joined %in% edge.names)]
        if (length(joined) > 0) 
            out[i, joined] <- 1
    }
    invisible(out)
}
fillPClust <-
function (cluster, k = "all", median.use = F, e.val.cutoff = NA, 
    p.val.cutoff = NA) 
{
    if (is.na(p.val.cutoff)) 
        p.val.cutoff <- get.iter.based.params("p.val.cutoff", 
            iter, params)
    if (is.na(e.val.cutoff)) 
        e.val.cutoff <- get.iter.based.params("e.val.cutoff", 
            iter, params)
    if (k == "all") {
        for (i in 1:cluster$k) cluster[[i]] <- fillPClust(cluster[[i]], 
            k = i, e.val.cutoff = e.val.cutoff, p.val.cutoff = p.val.cutoff)
        invisible(cluster)
    }
    else if (is.null(cluster$k) || (k <= cluster$k && k >= 1)) {
        m.out <- cluster$motif.out
        if (!is.null(m.out) && !is.null(m.out$e.values) && !all(is.na(m.out$e.values)) && 
            any(m.out$e.values < e.val.cutoff)) {
            cluster$p.clust <- get.p.clust(cluster, median.use = median.use, 
                p.val.cutoff = p.val.cutoff)
            cluster$e.val <- min(cluster$motif.out$e.values)
        }
        else {
            cluster$p.clust <- NA
            cluster$e.val <- min(cluster$motif.out$e.values)
        }
        invisible(cluster)
    }
}
find.best.weights.for.split <-
function (rows, row.gains, bad.move.temp, all.names, w.range = c(0.05, 
    0.5)) 
{
    not.rows <- all.names[!all.names %in% rows]
    fit.func <- function(w, rows, row.gains, bad.move.temp) {
        resp <- fit.gains.brlr(rows, row.gains, all.names = all.names, 
            weights = c(0.1, w))
        names(resp) <- names(row.gains)
        resp[rows] <- 1 - resp[rows]
        sample.prob <- exp((resp - 1)/bad.move.temp)
        y <- (sum(sample.prob[rows]) - sum(sample.prob[not.rows]))^2
        cat.new("FINDING BEST WEIGHTS:", w, y, "\n")
        return(y)
    }
    pars <- optimize(fit.func, w.range, tol = 0.01, rows = rows, 
        row.gains = row.gains, bad.move.temp = bad.move.temp)
    return(pars$minimum)
}
fit.gains.brlr <-
function (rows, row.gains, all.names = get.global("gene.ids"), 
    weights = c(1, 1), glm.safe = F) 
{
    gain.thresh <- 10
    w <- weights
    rm(weights)
    require(brlr)
    if (is.null(w)) {
        frac.1 <- length(rows)/length(all.names)
        frac.2 <- length(not.rows)/length(all.names)
    }
    else {
        frac.1 <- w[1]
        frac.2 <- w[2]
    }
    membership <- as.integer(all.names %in% rows)
    w <- membership * frac.2
    w[w == 0] <- frac.1
    names(membership) <- names(w) <- all.names
    if (glm.safe) {
        cat.new("over gains:\n")
        print.new(row.gains[row.gains > gain.thresh])
        row.gains[row.gains > gain.thresh] <- gain.thresh
    }
    fm <- brlr(membership ~ row.gains, weights = w)
    resp <- predict(fm, type = "response")
    return(resp)
}
gen.clust <-
function (rowNames, colNames, ratios, n.motifs = 2, fill = F) 
{
    c.tmp <- list()
    rowNames <- rowNames[which(rowNames %in% rownames(ratios))]
    colNames <- colNames[which(colNames %in% colnames(ratios))]
    c.tmp$nrows <- length(rowNames)
    c.tmp$ncols <- length(colNames)
    c.tmp$rows <- rowNames
    c.tmp$cols <- colNames
    if (fill) {
        c.tmp$k <- 999
        c.tmp$resid <- residOneClust(ratios[c.tmp$rows, c.tmp$cols], 
            varNorm = TRUE)
        c.tmp$motif.out <- motif.one.cluster(c.tmp)
        c.tmp$motif.out$mast.out <- NULL
        c.tmp <- fillPClust(c.tmp, k = 999, e.val.cutoff = e.val.cutoff, 
            p.val.cutoff = p.val.cutoff)
    }
    else {
        c.tmp$p.clust <- 1
        c.tmp$e.val <- rep(100, n.motifs)
    }
    return(c.tmp)
}
gen.clust.semirandom <-
function (ratios, clust.size = 15, ratios.thresh = 0.5, seed.genes = NULL, 
    seed.row = NULL) 
{
    if (is.null(seed.genes)) 
        seed.genes <- rownames(ratios)
    genes <- rownames(ratios)
    if (!is.null(seed.genes)) 
        genes <- seed.genes
    conds <- colnames(ratios)
    if (is.null(seed.row)) 
        seed.row <- sample(genes, 1)
    rats <- ratios[seed.row, ]
    while (!any(abs(rats) > ratios.thresh)) {
        seed.row <- sample(genes, 1)
        rats <- ratios[seed.row, ]
    }
    names(rats) <- conds
    cat.new("Seeding cluster using gene", seed.row, "and", length(seed.genes), 
        "genes\n")
    good.conds <- names(which(abs(rats) > ratios.thresh))
    cors <- apply(ratios[seed.genes, ], 1, cor, rats, use = "pair")
    best.rows <- genes[sort(cors, index = T, dec = T)$ix][1:clust.size]
    c.tmp <- gen.clust(best.rows, good.conds, ratios, fill = F)
    return(c.tmp)
}
gene.network.pvals <-
function (genes.test, cluster, networks = get.global("networks"), 
    net.q0 = get.global("net.q0"), net.names = names(net.q0), 
    gene.list = cluster$rows, force = F, is.cols = FALSE, min.edges = 1) 
{
    if (is.null(networks$adj.lists)) {
        networks$adj.lists <- list()
        for (i in net.names) networks$adj.lists[[i]] <- edge.list.to.adj.list(networks[[i]])
    }
    out <- list()
    for (i in net.names) {
        if ((!is.cols && i == "cond.sims") || (is.cols && i != 
            "cond.sims")) 
            next
        if (!force && (net.q0[i] <= 0 || is.null(networks[[i]]))) 
            next
        network <- networks[[i]]
        adj.list <- networks$adj.lists[[i]]
        n.in.out <- 0
        tmp <- length(unlist(adj.list[gene.list]))
        if (!is.na(tmp) && length(tmp) > 0) 
            n.in.out <- tmp
        n.in.in <- 0
        tmp <- unlist(adj.list[gene.list])
        if (length(tmp > 0)) 
            n.in.in <- sum(tmp %in% gene.list)
        out[[i]] <- rep(1, length(genes.test))
        names(out[[i]]) <- genes.test
        if (n.in.in + n.in.out <= 0) 
            next
        for (gene.test in genes.test) {
            tmp <- adj.list[[gene.test]]
            n.i.out <- length(tmp)
            if (n.i.out <= 0) 
                next
            n.i.in <- 0
            if (length(tmp) > 0) 
                n.i.in <- sum(gene.list %in% tmp)
            if (n.i.in <= 0) 
                next
            p.new <- phyper(n.i.in, n.i.in + n.in.in, n.i.out + 
                n.in.out, n.i.in + n.i.out, lower.tail = n.i.in >= 
                n.i.out)
            out[[i]][gene.test] <- p.new
        }
    }
    out
}
genes.in.how.many.clusters <-
function (clusters, genes = get.global("gene.ids"), up.to = clusters$k, 
    use.cols = F) 
{
    out <- integer(length(genes))
    names(out) <- genes
    g <- t(genes)
    if (!use.cols) {
        for (i in 1:up.to) out <- out + (g %in% clusters[[i]]$rows)
    }
    else {
        for (i in 1:up.to) out <- out + (g %in% clusters[[i]]$cols)
    }
    return(out)
}
genes.in.which.clusters <-
function (clusters, genes = get.global("gene.ids"), min = 1) 
{
    gene.names <- gene.coords$gene.name
    nms <- names(gene.names)
    out <- list()
    if (length(genes) <= 0) 
        return(out)
    for (i in 1:clusters$k) {
        glist <- character()
        for (g in genes) {
            if (!(g %in% nms) && g %in% gene.names) 
                g <- nms[which(gene.names == g)]
            grp <- grep(g, clusters[[i]]$rows, ignore.case = T)
            if (length(grp) > 0) 
                glist <- c(glist, clusters[[i]]$rows[grp])
        }
        if (length(glist) >= min) 
            out[[as.character(i)]] <- glist
    }
    return(out)
}
genes.not.in.any.clusters <-
function (clusterStack, all.genes = get.global("gene.ids")) 
{
    g <- character()
    for (i in 1:clusterStack$k) g <- c(g, clusterStack[[i]]$rows)
    g <- unique(g)
    out <- all.genes[!all.genes %in% g]
    return(out)
}
get.all.iter.based.params <-
function (iter = 1, params = get("params", .GlobalEnv), name = NULL, 
    verbose = F) 
{
    if (is.null(name)) {
        out <- list()
        out$iter <- iter
        for (name in names(params)) {
            if (substring(name, 1, 5) != "iter.") 
                next
            val <- get.iter.based.params(name, iter, params)
            real.name <- substring(name, 6)
            out[[real.name]] <- val
        }
        out$r0 <- out$resid.scaling
        out$p0 <- out$motif.scaling
        out$v0 <- out$volume.scaling
        out$q0 <- out$net.scaling
        out$net.q0 <- get.net.q0(params$net.max.weights, iter, 
            out$q0, verbose = verbose)
        return(out)
    }
    else {
        return(get.iter.based.params(name, iter, params))
    }
}
getBiMotifForCluster <-
function (cluster, gene.ids, n.motifs, min.seqs, max.seqs, min.width, 
    max.width, do.pal, seq.length, use.mcast = F, p.val.cutoff, 
    e.val.cutoff, motif.e.val.cutoff, meme.model = "zoops", mast.return = F, 
    unlink = T, verbose = F) 
{
    k <- cluster$k
    out <- list()
    n.genes <- length(gene.ids)
    p.val.mat <- rep(2, n.genes)
    names(p.val.mat) <- gene.ids
    out$p.values <- p.val.mat
    out$is.pal <- rep(do.pal, n.motifs)
    e.values <- out$e.values <- rep(e.val.cutoff, n.motifs)
    pssms <- out$pssms <- list()
    diagram.mat <- character(n.genes)
    names(diagram.mat) <- gene.ids
    out$too.few <- FALSE
    if (cluster$nrows < min.seqs) {
        if (verbose) 
            cat.new("Too few genes in clust ", k, ":", cluster$nrows, 
                "\n\n")
        out$too.few <- TRUE
        return(out)
    }
    if (verbose) 
        cat.new("Running MEME/MAST for cluster: ", k, "--------- --------- ------ ---- --- -- - - - -\n")
    genes <- cluster$rows
    sgenes <- genes[genes %in% gene.ids]
    if (is.null(cluster$opUpstream)) 
        seqs.500 <- opUpstream[sgenes]
    else seqs.500 <- cluster$opUpstream
    seqs.500 <- seqs.500[!is.na(seqs.500) & seqs.500 != ""]
    dups <- get.cluster.dup.seqs(cluster)
    sgenes <- sgenes[!dups]
    seqs.500 <- seqs.500[!dups]
    out$num.unique <- length(sgenes)
    if (length(sgenes) < min.seqs) {
        if (verbose) 
            cat.new("Too few unique genes in clust: ", k, "\n\n")
        out$too.few <- TRUE
        return(out)
    }
    else if (length(sgenes) > max.seqs) {
        if (verbose) 
            cat.new("Too many unique genes in clust: ", k, "\n\n")
        out$too.many <- TRUE
        return(out)
    }
    bgfname <- paste(tmp.prefix, sprintf("%03d", as.integer(k)), 
        ".bg", sep = "")
    bgseqs <- NULL
    if (is.null(upstream.bg)) {
        bgseqs <- opUpstream[!(gene.ids %in% genes)]
        bgdups <- duplicated(bgseqs)
        bgseqs <- bgseqs[!bgdups]
        bgseqs <- bgseqs[!(bgseqs %in% seqs.500)]
    }
    else {
        bad.res <- c("Y", "R", "K", "W", "N", "S")
        whch <- unlist(sapply(bad.res, function(i) grep(i, names(upstream.bg))), 
            use.names = F)
        if (length(whch) > 0) 
            upstream.bg <- upstream.bg[-whch]
    }
    seqs <- seqs.500
    seqs <- seqs[!is.na(seqs) & seqs != ""]
    if (seq.length > 0 && seq.length < nchar(seqs[1])) 
        seqs <- subSeq(seqs.500, len = seq.length, direction = "back")
    fname <- paste(tmp.prefix, sprintf("%03d", as.integer(k)), 
        ".fst", sep = "")
    meme.out <- list()
    meme.info <- memeP.info <- NULL
    e.cutoff <- e.val.cutoff
    if (do.pal) {
        meme.out <- runMemePalNonPal(sgenes, seqs, bgseqs = bgseqs, 
            fname = fname, bgfname = bgfname, nmotif = n.motifs, 
            e.value.cutoff = e.cutoff, bg.list = upstream.bg, 
            minw = min.width, maxw = max.width, model = meme.model, 
            unlink = unlink, verbose = verbose)
        meme.info <- getMemeMotifInfo(meme.out$nonpal)
        tries <- 1
        while (tries < 20 && length(meme.info) < n.motifs) {
            meme.info <- getMemeMotifInfo(meme.out$nonpal)
            tries <- tries + 1
        }
        for (mot in 1:length(meme.info)) if (verbose) 
            cat.new("E-val non pal: (", mot, ")", meme.info[[mot]]$e.value, 
                ";\tnum seqs", length(sgenes), "\n")
        memeP.info <- getMemeMotifInfo(meme.out$pal)
        tries <- 1
        while (tries < 20 && length(memeP.info) < n.motifs) {
            memeP.info <- getMemeMotifInfo(meme.out$pal)
            tries <- tries + 1
        }
        for (mot in 1:length(meme.info)) if (verbose) 
            cat.new("E-val pal: (", mot, ")", memeP.info[[mot]]$e.value, 
                ";\tnum seqs", length(sgenes), "\n")
    }
    else {
        meme.out$nonpal <- runMeme(sgenes, seqs, bgseqs = bgseqs, 
            fname = fname, bgfname = bgfname, nmotif = n.motifs, 
            e.value.cutoff = e.cutoff, bg.list = upstream.bg, 
            minw = min.width, maxw = max.width, model = meme.model, 
            unlink = unlink, verbose = verbose)
        meme.info <- getMemeMotifInfo(meme.out$nonpal)
        tries <- 1
        while (tries < 20 && length(meme.info) < n.motifs) {
            meme.info <- getMemeMotifInfo(meme.out$nonpal)
            tries <- tries + 1
        }
        if (length(meme.info) > 0) {
            for (mot in 1:length(meme.info)) if (verbose) 
                cat.new("E-val (non pal only): (", mot, ")", 
                  meme.info[[mot]]$e.value, ";\tnum seqs", length(sgenes), 
                  "\n")
        }
        else {
            if (verbose) 
                cat.new("E-val (non pal only): ", e.cutoff, "; num genes ", 
                  cluster$nrows, "\n")
        }
    }
    meme.todo <- list(meme.out$nonpal)
    if (length(meme.info) > 0) {
        if (verbose) 
            cat.new("Got", length(meme.info), "good motifs for cluster", 
                k, "\n")
        pssms <- list()
        for (i in 1:length(meme.info)) {
            e.values[i] <- e.cutoff
            pssms[[i]] <- NULL
            if (meme.info[[i]]$e.value < e.cutoff || (!is.null(memeP.info) && 
                memeP.info[[i]]$e.value < e.cutoff)) {
                if (do.pal) {
                  if (meme.info[[i]]$e.value <= memeP.info[[i]]$e.value) {
                    if (verbose) 
                      cat.new("Using non--palindromic motif for motif", 
                        i, ", cluster: ", k, "\n")
                    e.values[i] <- meme.info[[i]]$e.value
                    pssms[[i]] <- meme.info[[i]]$pssm
                  }
                  else {
                    if (verbose) 
                      cat.new("Using palindromic motif for motif", 
                        i, ", cluster: ", k, "\n")
                    meme.todo[[2]] <- meme.out$pal
                    out$is.pal[i] <- TRUE
                    e.values[i] <- memeP.info[[i]]$e.value
                    pssms[[i]] <- memeP.info[[i]]$pssm
                  }
                }
                else {
                  e.values[i] <- meme.info[[i]]$e.value
                  pssms[[i]] <- meme.info[[i]]$pssm
                }
            }
        }
        e.cutoff <- e.val.cutoff
        motif.e.cutoff <- motif.e.val.cutoff
        if (!any(e.values <= motif.e.cutoff)) 
            motif.e.cutoff <- min(e.values) + 1
        if (verbose) 
            cat.new("Using", sum(e.values <= motif.e.cutoff), 
                "motifs for MAST/MCAST search with E-value <=", 
                motif.e.cutoff, "\n")
        mast.fname <- paste(tmp.prefix, "_mast_", sprintf("%03d", 
            as.integer(k)), ".fst", sep = "")
        mast.bgfname <- paste(tmp.prefix, "_mast_", sprintf("%03d", 
            as.integer(k)), ".bg", sep = "")
        memeOutFname <- paste(tmp.prefix, sprintf("%03d", as.integer(k)), 
            ".meme.out", sep = "")
        mast.out <- mcast.out <- list()
        search.seqs <- opUpstream[gene.ids]
        search.seqs <- search.seqs[!is.na(search.seqs) & search.seqs != 
            ""]
        if (!use.mcast) {
            for (ind in length(meme.todo)) {
                mast.out[[ind]] <- runMast(meme.todo[[ind]], 
                  names(search.seqs), search.seqs, bgseqs = bgseqs, 
                  bg.list = upstream.bg, fname = mast.fname, 
                  memeOutFname = memeOutFname, bgfname = mast.bgfname, 
                  p.value.cutoff = p.val.cutoff, e.value.cutoff = e.cutoff, 
                  motif.e.value.cutoff = motif.e.cutoff, verbose = verbose, 
                  unlink = unlink)
            }
            pval.eval <- getMastPValuesAndEValues(mast.out[[1]], 
                genes)
            out$mast.info <- list()
            for (g in genes) out$mast.info[[g]] <- pval.eval[[g]]
            p.vals <- unlist(lapply(pval.eval, function(i) i$p.value))
            genes.passed <- names(p.vals)
            p.val.mat[genes.passed] <- p.vals[genes.passed]
            diagrams <- unlist(lapply(pval.eval, function(i) i$diagram))
            diagram.mat[genes.passed] <- diagrams[genes.passed]
        }
        else {
            p.v.cutoff <- p.val.cutoff
            if (p.v.cutoff > 0.01) 
                p.v.cutoff <- 0.01
            for (ind in length(meme.todo)) {
                mcast.out[[ind]] <- runMCast(meme.todo[[ind]], 
                  gene.ids, search.seqs, bgseqs = bgseqs, bg.list = upstream.bg, 
                  fname = mast.fname, memeOutFname = memeOutFname, 
                  bgfname = mast.bgfname, p.value.cutoff = p.v.cutoff, 
                  e.value.cutoff = e.cutoff, max.gap = max.motif.gap, 
                  motif.e.value.cutoff = motif.e.cutoff, verbose = verbose, 
                  unlink = unlink)
            }
            mcast.info <- parse.mcast.output(mcast.out[[1]], 
                parse = genes, verbose = F, p.val.cutoff = p.v.cutoff)
            if (!is.null(mcast.info)) {
                p.vals <- mcast.info$hits[["p-value"]]
                p.vals[p.vals <= 0] <- 1e-12
                names(p.vals) <- mcast.info$hits$ID
                genes.passed <- names(p.vals)
                p.val.mat[genes.passed] <- p.vals[genes.passed]
                tmp.seqs <- opUpstream[genes]
                diagrams <- mcast.info.to.diagrams(mcast.info$info, 
                  meme.info = meme.info, seqs = tmp.seqs)
                if (!is.null(diagrams) && length(diagrams) > 
                  0) 
                  diagram.mat[genes] <- diagrams[genes]
                out$mast.info <- list()
                for (g in genes) out$mast.info[[g]] <- mcast.info$info[[g]]
            }
        }
        out$p.values <- log10(p.val.mat)
        out$diagrams <- diagram.mat
        out$e.values <- e.values
        out$pssms <- pssms
        out$meme.out <- meme.out
        if (mast.return) {
            if (!use.mcast) 
                out$mast.out <- mast.out
            else out$mast.out <- mcast.out
        }
        if (use.mcast) 
            out$use.mcast <- TRUE
    }
    else {
        if (verbose) 
            cat.new("No motif by any method had a good enough E-val for cluster ", 
                k, "\n")
    }
    gc()
    invisible(out)
}
get.cluster.dup.seqs <-
function (cluster) 
{
    rows <- cluster$rows
    seqs <- get.cluster.seqs(cluster)
    is.dup.seq <- duplicated(seqs[rows])
    names(is.dup.seq) <- rows
    return(is.dup.seq)
}
get.cluster.dup.seqs.GOOD <-
function (cluster) 
{
    rows <- cluster$rows
    homs <- dup.seqs[dup.seqs[, 1] %in% rows, ]
    homs2 <- homs[homs[, 2] %in% rows, ]
    if (nrow(homs2) > 0) {
        for (i in 1:nrow(homs2)) rows[rows == homs2[i, 2]] <- homs2[i, 
            1]
    }
    is.dup.seq <- duplicated(rows)
    names(is.dup.seq) <- cluster$rows
    return(is.dup.seq)
}
get.cluster.likelihood <-
function (cluster, resp) 
{
    score <- sum(log10(1 - resp[cluster$rows]) + log10(resp[gene.ids[!(gene.ids %in% 
        cluster$rows)]]))
    return(score)
}
get.cluster.seqs <-
function (cluster) 
{
    rows <- cluster$rows
    seqs <- opUpstream[rows]
    return(seqs)
}
get.cluster.subnet.p.values <-
function (cluster, networks, gene.ids, edge.centric = T, net.q0 = NA, 
    force = F, rows.only = F, cols.only = F) 
{
    out <- numeric()
    if (!cols.only) 
        out <- subnet.p.values(cluster$rows, networks, gene.ids, 
            edge = edge.centric, net.q0, force)
    if (!rows.only) {
        if (force || (!is.na(net.q0) && !is.na(net.q0["cond.sims"]) && 
            net.q0["cond.sims"] > 0)) 
            out["cond.sims"] <- condition.sim.pvalue(cluster$cols, 
                networks, all.conditions, edge = edge.centric)
        else out["cond.sims"] <- 2
    }
    else out["cond.sims"] <- 2
    out
}
get.col.gains.for.cluster <-
function (cluster, all.conds = get.global("all.conditions"), 
    ratios = get.global("ratios"), networks = get.global("networks"), 
    net.q0 = get.global("net.q0"), response = T, bad.move.temp = get.global("bad.move.temperature"), 
    r.sig = get.global("r.sig"), plot.it = F) 
{
    cols <- cluster$cols
    vars <- get.vars.for.cluster(cluster, opt = "cols", genes = all.conds, 
        ratios = ratios, var.norm = T, r.sig = r.sig)
    if (q0 > 0 && net.q0["cond.sims"] > 0) {
        net.pvals <- get.network.pvals.for.cluster(cluster, net.names = "cond.sims", 
            gene.ids = all.conds, networks = networks, net.q0 = net.q0, 
            ignore = NULL)
        if (sum(net.pvals) > 0) {
            net.pvals <- net.pvals - mean(net.pvals[, cols], 
                na.rm = T)
            tmp.sd <- sd(net.pvals[, cols], na.rm = T)
            if (!is.na(tmp.sd) && tmp.sd != 0) 
                net.pvals <- net.pvals/tmp.sd
            net.pvals <- net.pvals * net.q0["cond.sims"]
        }
        gains <- r0 * vars + q0 * as.numeric(net.pvals)
    }
    else {
        gains <- r0 * vars
    }
    tmp.sd <- sd(gains[cols], na.rm = T)
    if (!is.na(tmp.sd) && tmp.sd != 0) 
        gains <- gains/tmp.sd
    gains <- gains - mean(gains[cols], na.rm = T)
    gains <- -gains
    if (response) {
        resp <- pnorm(gains, mean = 1)
        names(resp) <- names(gains)
        resp[cols] <- 1 - pnorm(gains[cols], mean = -3)
        if (plot.it) {
            hist(gains[cols], breaks = 20, main = "gains[cols]")
            q.hist <- hist(gains, breaks = 20, main = "col.gains")
            points(gains, resp * max(q.hist$counts) * 0.9)
            points(gains[cols], resp[cols] * max(q.hist$counts) * 
                0.9, col = "red")
        }
        return(resp)
    }
    return(gains)
}
get.cors.for.cluster <-
function (cluster, genes = get.global("gene.ids"), ratios = get.global("ratios")) 
{
    rat <- ratios[cluster$rows, cluster$cols, drop = F]
    prof <- apply(rat, 2, mean, na.rm = T)
    cors <- apply(ratios[genes, cluster$cols, drop = F], 1, cor, 
        prof, use = "pairwise")
    names(cors) <- genes
    return(cors)
}
getEntropy <-
function (pssm) 
{
    win.size <- dim(pssm)[1]
    pssm[pssm == 0] <- 1e-05
    entropy <- vector(mode = "numeric", length = win.size)
    for (i in 1:win.size) {
        entropy[i] <- -1 * sum(pssm[i, ] * log2(pssm[i, ]))
    }
    return(entropy)
}
get.global <-
function (name) 
{
    return(get(name, .GlobalEnv))
}
get.iter.based.params <-
function (name, iter, params = get("params", .GlobalEnv)) 
{
    if (is.null(params[[name]])) 
        name <- paste("iter.", name, sep = "")
    out <- params[[name]][1]
    if (length(params[[name]]) > 1) {
        tmp <- vector()
        tmp[as.integer(names(params[[name]]))] <- params[[name]]
        out <- last.element(tmp[1:iter][!is.na(tmp[1:iter])])
    }
    return(out)
}
get.log <-
function (prefix = "biclust.log.%03d.txt", output.dir = "output/") 
{
    if (exists("out.log")) 
        return(out.log)
    if (!file.exists(output.dir)) 
        try(dir.create(output.dir, recursive = T, showWarnings = F))
    for (i in 1:999) {
        f <- paste(output.dir, sprintf(prefix, as.integer(i)), 
            sep = "")
        if (file.exists(f)) 
            next
        f <- paste(getwd(), f, sep = "/")
        out.log <<- f
        cat.new("Logging all output to", out.log, "\n")
        break
    }
    return(out.log)
}
get.loglik <-
function (rows, row.gains, mean = -1, sd = 1) 
{
    not.rows <- gene.ids[!(gene.ids %in% rows)]
    resp <- pnorm(row.gains, mean = mean, sd = sd)
    names(resp) <- names(row.gains)
    frac.1 <- length(rows)/length(gene.ids)
    frac.2 <- length(not.rows)/length(gene.ids)
    return(sum(log10(1 - resp[rows]))/frac.1 + sum(log10(resp[not.rows]))/frac.2)
}
get.mast.pvals <-
function (mast.output, in.genes = NULL) 
{
    start <- grep("SECTION III: ANNOTATED SEQUENCES", mast.output)
    end <- grep("\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*\\*", 
        mast.output[(start + 3):length(mast.output)]) + start + 
        3
    line.starts <- grep("LENGTH = ", mast.output[(start + 2):(start + 
        1 + end)]) + start + 1
    residues <- c("G", "A", "T", "C")
    out <- list()
    for (i in 1:length(line.starts)) {
        l <- line.starts[i]
        gene <- mast.output[l - 2]
        if (!is.null(in.genes) && !gene %in% in.genes) 
            next
        out[[gene]] <- list()
        i.tmp <- 2
        l.first = ""
        while (l.first == "") {
            i.tmp <- i.tmp + 1
            l.first <- mast.output[l + i.tmp]
        }
        additional <- i.tmp - 3
        if (i >= length(line.starts)) 
            l.next <- end
        else l.next <- line.starts[i + 1] - 2
        l.last <- mast.output[l.next - 3]
        if (l.next - l <= 5) {
            next
        }
        seq.inds <- seq(l + 7 + additional, l.next - 3, by = 6)
        pos.offset <- as.integer(strsplit(mast.output[seq.inds[1]], 
            " ")[[1]][1])
        first.nonspace.ind <- which(strsplit(mast.output[seq.inds[1]], 
            "")[[1]] %in% residues)[1]
        pos.offsets <- as.integer(substring(mast.output[seq.inds], 
            1, first.nonspace.ind - 1))
        motif.inds <- seq(l + 3 + additional, l.next - 3, by = 6)
        lines.tmp <- space.pad(mast.output[motif.inds], 80)
        mots <- paste(substring(lines.tmp, first.nonspace.ind), 
            sep = "", collapse = "")
        mots <- gsub("\\[", " ", mots)
        mots <- gsub("\\]", " ", mots)
        mots <- gsub("<", " ", mots)
        mots <- gsub(">", " ", mots)
        mots <- strsplit(mots, "\\s+", perl = T)[[1]]
        mots <- mots[mots != ""]
        mots <- as.integer(mots)
        pval.inds <- seq(l + 4 + additional, l.next - 3, by = 6)
        lines.tmp <- space.pad(mast.output[pval.inds], 80)
        tmp <- paste(substring(lines.tmp, first.nonspace.ind), 
            sep = "", collapse = "")
        pvals <- as.numeric(strsplit(tmp, "\\s+", perl = T)[[1]])
        pvals <- pvals[!is.na(pvals)]
        tmp <- gsub("[^\\s]", "X", tmp, perl = T)
        posns <- grep("\\sX", substring(tmp, 1:(nchar(tmp) - 
            1), 2:nchar(tmp)), perl = T) + 1
        posns <- posns + pos.offset - 1
        out[[gene]]$pvals <- pvals
        out[[gene]]$mots <- mots
        out[[gene]]$posns <- posns
    }
    return(out)
}
getMastPValuesAndEValues <-
function (mastOutput, get.p.values = NULL) 
{
    out <- list()
    lines <- grep("COMBINED P-VALUE", mastOutput)
    if (length(lines) > 0) {
        for (i in 1:length(lines)) {
            gene <- mastOutput[lines[i] - 2]
            splitted <- unique(unlist(strsplit(mastOutput[lines[i]], 
                " ")))
            p.val <- as.numeric(splitted[7])
            e.val <- as.numeric(splitted[9])
            splitted <- unique(unlist(strsplit(mastOutput[lines[i] + 
                1], " ")))
            diagram <- splitted[3]
            j <- 2
            while (substring(diagram, nchar(diagram), nchar(diagram)) == 
                "_") {
                splitted <- unique(unlist(strsplit(mastOutput[lines[i] + 
                  j], " ")))
                diagram <- paste(diagram, splitted[2], sep = "")
                j <- j + 1
            }
            out[[gene]] <- list(p.value = p.val, e.value = e.val, 
                diagram = diagram)
        }
    }
    else {
        out[["none"]] <- list(p.value = NA, e.value = NA, diagram = character(0))
    }
    if (!is.null(get.p.values)) {
        tmp <- get.mast.pvals(mastOutput, get.p.values)
        for (g in names(tmp)) {
            out[[g]]$pvals <- tmp[[g]]$pvals
            out[[g]]$posns <- tmp[[g]]$posns
            out[[g]]$mots <- tmp[[g]]$mots
        }
    }
    return(out)
}
getMemeMotifInfo <-
function (memeOutput) 
{
    out <- list()
    lines <- grep("^MOTIF\\s+\\d", memeOutput, perl = T)
    if (length(lines) <= 0) 
        lines <- grep("^MOTIF\\s+", memeOutput, perl = T)
    if (length(lines) > 0) {
        pssms <- getMemeMotifPssm(memeOutput, n.motif = length(lines))
        for (i in 1:length(lines)) {
            line <- memeOutput[lines[i]]
            splitted <- unlist(strsplit(line, "[\\t\\s]", perl = T))
            splitted <- splitted[splitted != ""]
            motif <- as.integer(splitted[2])
            width <- as.integer(splitted[5])
            sites <- as.integer(splitted[8])
            llr <- as.integer(splitted[11])
            e.value <- as.numeric(sub("\\+", "", splitted[14]))
            pssm <- pssms[[motif]]$pssm
            posns <- list()
            l2 <- grep(paste("Motif", motif, "sites sorted by position p-value"), 
                memeOutput) + 4
            while (memeOutput[l2] != "--------------------------------------------------------------------------------") {
                line <- unlist(strsplit(memeOutput[l2], "[\\t\\s]", 
                  perl = T))
                line <- line[line != ""]
                posns[[line[1]]] <- list(strand = line[2], start = line[3], 
                  p.value = line[4], site = line[6])
                l2 <- l2 + 1
            }
            out[[motif]] <- list(width = width, sites = sites, 
                llr = llr, e.value = e.value, pssm = pssm, posns = posns)
        }
    }
    out
}
getMemeMotifPssm <-
function (memeOut, n.motif = 1) 
{
    pssms <- list()
    for (i in 1:n.motif) {
        m.start <- paste("Motif ", i, " position-specific probability matrix", 
            sep = "")
        m.line1 <- grep(m.start, memeOut)
        if (length(m.line1) > 0) {
            m.desc <- unlist(strsplit(memeOut[m.line1 + 2], " "))
            winLen <- as.numeric(m.desc[6])
            e.val <- as.numeric(m.desc[10])
            pssm <- matrix(, nrow = winLen, ncol = 4)
            for (j in 1:winLen) {
                pssm[j, ] <- as.numeric(unlist(strsplit(memeOut[m.line1 + 
                  2 + j], " "))[c(2, 4, 6, 8)])
            }
            tmp.p <- list()
            tmp.p$pssm <- pssm
            tmp.p$e.val <- e.val
            pssms[[i]] <- tmp.p
            rm(tmp.p)
        }
        else {
            tmp.p <- list()
            tmp.p$pssm <- NULL
            tmp.p$e.val <- 99999
            pssms[[i]] <- tmp.p
        }
    }
    return(pssms)
}
get.motif.pvals.for.cluster <-
function (cluster, genes = get.global("gene.ids"), var.norm = T, 
    p.sig = 0.1) 
{
    if (!is.null(cluster$motif.out) && !is.null(cluster$motif.out$p.values)) {
        pv <- cluster$motif.out$p.values
        pvals <- pv[genes]
    }
    else {
        pvals <- rep(0, length(genes))
    }
    names(pvals) <- genes
    if (var.norm) {
        pvals <- pvals - mean(pvals, na.rm = T)
        clust.var <- sd(pvals[cluster$rows], na.rm = T)
        if (!is.na(clust.var) && clust.var != 0) 
            pvals <- pvals/(clust.var + p.sig)
    }
    return(pvals)
}
get.net.q0 <-
function (net.max.weights, iter, in.net.q0, verbose = F) 
{
    out.q0 <- net.max.weights
    for (i in names(net.max.weights)) {
        if (verbose) 
            cat.new("Network weight", i, ":", out.q0[i], "\n")
        out.q0[i] <- out.q0[i] * in.net.q0
        if (verbose) 
            cat.new("Adjusted network mixing weight", i, ":", 
                out.q0[i], "\n")
    }
    return(out.q0)
}
get.network.pvals.for.cluster <-
function (cluster, net.names = "all", gene.ids = get.global("gene.ids"), 
    networks = get.global("networks"), net.q0 = get.global("net.q0"), 
    ignore = c("cond.sims", "COG.code")) 
{
    if (net.names == "all") 
        net.names <- names(net.q0)
    if (!is.null(ignore)) 
        net.names <- net.names[!(net.names %in% ignore)]
    net.names <- net.names[!(net.names %in% names(which(net.q0 == 
        0)))]
    is.cols <- FALSE
    if (net.names == "cond.sims" && !"cond.sims" %in% ignore) 
        is.cols <- TRUE
    tmp <- gene.network.pvals(gene.ids, cluster, networks, net.q0, 
        net.names = net.names, force = F, is.cols = is.cols)
    if (length(tmp[[1]]) > 1) {
        net.pvals <- matrix(unlist(tmp), nrow = length(tmp[[1]]))
        rownames(net.pvals) <- names(tmp[[1]])
        colnames(net.pvals) <- names(tmp)
    }
    else if (length(tmp[[1]]) == 1) {
        net.pvals <- unlist(tmp)
        names(net.pvals) <- gene.ids
    }
    else if (length(tmp[[1]]) == 0) {
        net.pvals <- numeric(length(gene.ids))
        names(net.pvals) <- gene.ids
    }
    return(t(net.pvals))
}
get.network.stats <-
function (clusterStack, networks = get("networks", .GlobalEnv), 
    verbose = F, mean.type = 2, net.q0 = NA) 
{
    out <- rep(2, length(networks) - 4)
    names <- character()
    if (is.null(clusterStack[[1]]$net.p.old)) 
        return(out)
    vals <- t(sapply(clusterStack[1:clusterStack$k], "[[", "net.p.old"))
    vals <- log10(vals)
    names(out) <- colnames(vals)
    if (!verbose) {
        if (mean.type != 2) 
            vals[vals < -98 | vals >= log10(2)] <- NA
        else vals[vals < -98 | vals > log10(2)] <- NA
        if (all(is.na(vals))) 
            return(out)
        mean.q <- apply(vals, 2, mean, na.rm = T)
        mean.q[is.na(mean.q)] <- 1
        return(mean.q)
    }
    else {
        vv1 <- vals
        vv1[vv1 < -98 | vv1 >= log10(2)] <- NA
        var.q <- apply(vv1, 2, var, na.rm = T)
        var.q[is.na(var.q)] <- 1
        var.q[var.q == 0] <- 1
        mean.q <- apply(vv1, 2, mean, na.rm = T)
        mean.q[is.na(mean.q)] <- 1
        vv2 <- vals
        vv2[vv2 < -98 | vv2 > log10(2)] <- NA
        var.q2 <- apply(vv2, 2, var, na.rm = T)
        var.q2[is.na(var.q2)] <- 1
        var.q2[var.q2 == 0] <- 1
        mean.q2 <- apply(vv2, 2, mean, na.rm = T)
        mean.q2[is.na(mean.q2)] <- 1
        for (i in names(mean.q)) {
            cat.new("Network", i, ": mean =", mean.q[i], "+/-", 
                var.q[i], "\t\tmean2 =", mean.q2[i], "+/-", var.q2[i], 
                "\n")
        }
        if (mean.type == 2) 
            return(mean.q2)
        else return(mean.q)
    }
}
get.organism <-
function (choices = c("halo", "hpy", "yeast", "ecoli", "bsubt")) 
{
    choice <- menu(choices, title = "Choose your organism:")
    if (choice == 0) 
        choice <- 1
    return(choices[choice])
}
get.param <-
function (name, params = get("params", .GlobalEnv)) 
{
    return(params[[name]])
}
get.p.clust <-
function (cluster, median.use = F, rows = cluster$rows, p.val.cutoff = NA) 
{
    if (is.na(p.val.cutoff)) 
        p.val.cutoff <- get.iter.based.params("p.val.cutoff", 
            iter, params)
    tmp.mat <- cluster$motif.out$p.values[rows]
    tmp.mat[tmp.mat >= log10(p.val.cutoff)] <- NA
    p.clust <- NA
    if (!median.use) 
        p.clust <- mean(tmp.mat, na.rm = T)
    else p.clust <- median(tmp.mat, na.rm = T)
    return(p.clust)
}
get.row.gains.for.cluster <-
function (cluster, gene.ids = get.global("gene.ids"), ratios = get.global("ratios"), 
    networks = get.global("networks"), net.q0 = get.global("net.q0"), 
    response = T, bad.move.temp = get.global("bad.move.temperature"), 
    r.sig = get.global("r.sig"), weight.fixed = 0.5, find.best.weights = F, 
    plot.it = F) 
{
    rows <- cluster$rows
    vars <- get.vars.for.cluster(cluster, opt = "rows", genes = gene.ids, 
        ratios = ratios, var.norm = T, r.sig = r.sig)
    pvals <- get.motif.pvals.for.cluster(cluster, genes = gene.ids, 
        var.norm = T, p.sig = 1)
    use.networks <- (q0 > 0 && !is.null(net.q0) && sum(net.q0) > 
        0)
    if (use.networks) {
        net.pvals <- get.network.pvals.for.cluster(cluster, net.names = "all", 
            gene.ids = gene.ids, networks = networks, net.q0 = net.q0, 
            ignore = c("cond.sims", "COG.code"))
        net.pvals <- t(log10(net.pvals + 0.001))
        net.means <- apply(net.pvals[rows, ], 2, mean, na.rm = T)
        net.pvals <- t(apply(net.pvals, 1, "-", net.means))
        net.sds <- apply(net.pvals[rows, ], 2, sd, na.rm = T)
        net.sds[net.sds == 0] <- 1
        q.sig <- 0.3
        net.pvals <- t(apply(net.pvals, 1, "/", (net.sds + q.sig)))
        tmp.q0 <- net.q0[colnames(net.pvals)]/sum(net.q0[colnames(net.pvals)])
        net.gains <- t(apply(net.pvals, 1, "*", tmp.q0))
        net.gains <- apply(net.gains, 1, sum)
        tmp.sd <- sd(net.gains[rows], na.rm = T)
        if (tmp.sd == 0) 
            tmp.sd <- 1
        net.gains <- net.gains/(tmp.sd + q.sig)
        gains <- r0 * vars + p0 * pvals + q0 * net.gains
    }
    else {
        gains <- r0 * vars + p0 * pvals
    }
    tmp.sd <- sd(gains[rows], na.rm = T)
    if (!is.na(tmp.sd) && tmp.sd != 0) 
        gains <- gains/tmp.sd
    gains <- gains - mean(gains[rows], na.rm = T)
    gains <- -gains
    if (response) {
        if (find.best.weights) {
            ww <- find.best.weights.for.split(rows, gains, bad.move.temp, 
                all.names = gene.ids, w.range = c(0.05, 0.75))
        }
        else {
            ww <- weight.fixed
        }
        if (plot.it) 
            cat.new("Using WEIGHT:", ww, "\n")
        resp <- fit.gains.brlr(rows, gains, all.names = gene.ids, 
            weights = c(ww, 1))
        names(resp) <- names(gains)
        resp[rows] <- 1 - resp[rows]
        if (plot.it) {
            plot(vars, pvals, main = paste("Cluster:", cluster$k))
            points(vars[rows], pvals[rows], col = "red", pch = 20)
            hist(gains[rows], breaks = 20, main = "gains[rows]")
            q.hist <- hist(gains, breaks = 50, main = "row.gains")
            points(gains, resp * max(q.hist$counts) * 0.9)
            points(gains[rows], resp[rows] * max(q.hist$counts) * 
                0.9, col = "red")
        }
        return(resp)
    }
    return(gains)
}
get.vars.for.cluster <-
function (cluster, genes = get.global("gene.ids"), opt = c("rows", 
    "cols"), ratios = get.global("ratios"), var.norm = T, r.sig = get.global("r.sig"), 
    allow.anticor = get.global("allow.anticor")) 
{
    opt <- match.arg(opt)
    rows <- cluster$rows
    cols <- cluster$cols
    if (opt == "rows") {
        r <- ratios[rows, cols]
        r.all <- ratios[genes, cols]
        avg.rows <- apply(r, 2, mean, na.rm = T)
        devs <- apply(r.all, 1, "-", avg.rows)
        vars <- apply(devs, 2, var, na.rm = T)
        vars <- log10(vars)
        if (allow.anticor) {
            devs.2 <- apply(r.all, 1, "-", -avg.rows)
            vars.2 <- apply(devs.2, 2, var, na.rm = T)
            vars.2 <- log10(vars.2)
            vars <- cbind(vars, vars.2)
            vars <- apply(vars, 1, min)
        }
        if (var.norm) {
            vars <- vars - mean(vars[rows], na.rm = T)
            tmp.sd <- sd(vars[rows], na.rm = T)
            if (!is.na(tmp.sd) && tmp.sd != 0) 
                vars <- vars/(tmp.sd + r.sig)
        }
        return(vars)
    }
    else {
        r.all <- ratios[rows, ]
        vars <- log10(apply(r.all, 2, var, na.rm = T)/abs(apply(r.all, 
            2, mean, na.rm = T)))
        names(vars) <- colnames(ratios)
        if (var.norm) {
            vars <- vars - mean(vars[cluster$cols], na.rm = T)
            tmp.sd <- sd(vars[cluster$cols], na.rm = T)
            if (!is.na(tmp.sd) && tmp.sd != 0) 
                vars <- vars/(tmp.sd + r.sig)
        }
        return(vars)
    }
}
get.vizmap <-
function (nodes, edges, nsize = NULL, nshape = NULL, ncol = NULL, 
    nfill = NULL, ntyp = NULL, lcol = NULL, lfont = NULL, lsize = NULL, 
    loff.x = NULL, loff.y = NULL, ecol = NULL, etyp = NULL, ewid = NULL, 
    arrow = NULL, asize = NULL, opt = c("default", "random")) 
{
    node.count <- length(nodes)
    edge.count <- 0
    if (!is.null(edges)) 
        edge.count <- nrow(edges)
    opt <- match.arg(opt)
    if (opt == "random") {
        if (is.null(nshape)) 
            nshape <- sample(21:25, node.count, rep = T)
        if (is.null(nsize)) 
            nsize <- runif(node.count) * 2 + 0.5
        if (is.null(ncol)) 
            ncol <- sample(node.count)
        if (is.null(nfill)) 
            nfill <- sample(node.count)
        if (is.null(ntyp)) 
            ntyp <- sample(node.count)
        if (is.null(lcol)) 
            lcol <- sample(node.count)
        if (is.null(lfont)) 
            lfont <- sample(node.count)
        if (is.null(lsize)) 
            lsize <- runif(node.count) * 2 + 0.5
        if (is.null(loff.x)) 
            loff.x <- 0
        if (is.null(loff.y)) 
            loff.y <- 0
        if (is.null(ecol)) 
            ecol <- sample(edge.count)
        if (is.null(etyp)) 
            etyp <- sample(edge.count)
        if (is.null(ewid)) 
            ewid <- (runif(edge.count) + 0.3) * 3
        if (is.null(arrow)) 
            arrow <- rep(2, edge.count)
        if (is.null(asize)) 
            asize <- (runif(edge.count) + 0.2) * 0.1
    }
    else if (opt == "default") {
        if (is.null(nshape)) 
            nshape <- 21
        if (is.null(nsize)) 
            nsize <- 5
        if (is.null(ncol)) 
            ncol <- "black"
        if (is.null(nfill)) 
            nfill <- "darkgrey"
        if (is.null(ntyp)) 
            ntyp <- 1
        if (is.null(lcol)) 
            lcol <- "white"
        if (is.null(lfont)) 
            lfont <- 1
        if (is.null(lsize)) 
            lsize <- 1
        if (is.null(loff.x)) 
            loff.x <- 0
        if (is.null(loff.y)) 
            loff.y <- 0
        if (is.null(ecol)) 
            ecol <- rep("darkgrey", edge.count)
        if (is.null(etyp)) 
            etyp <- rep(1, edge.count)
        if (is.null(ewid)) 
            ewid <- rep(1, edge.count)
        if (is.null(arrow)) 
            arrow <- rep(0, edge.count)
        if (is.null(asize)) 
            asize <- rep(0, edge.count)
    }
    out <- list(nodes = nodes, edges = edges, nsize = nsize, 
        nshape = nshape, ncol = ncol, nfill = nfill, lcol = lcol, 
        lfont = lfont, lsize = lsize, ntyp = ntyp, loff.x = loff.x, 
        loff.y = loff.y, ecol = ecol, etyp = etyp, ewid = ewid, 
        arrow = arrow, asize = asize)
    return(out)
}
graphics.to.postscript <-
function (fname = "Rplots.ps") 
{
    dev.copy(device = postscript, paper = "letter")
    dev.off()
    if (fname != "Rplots.ps") 
        system(paste("mv -fv Rplots.ps", fname))
    system(paste("gzip -fv ", fname))
    system(paste("gv ", fname, ".gz &", sep = ""))
}
histPClust <-
function (clust) 
{
    return(unlist(lapply(clust, function(i) i$p.clust), use.names = F))
}
histResid <-
function (clust, n1 = 0, n2 = 0) 
{
    return(unlist(lapply(clust, function(i) i$resid), use.names = F))
}
init <-
function (is.cluster.node = F) 
{
    initialize.biclust(T, is.cluster.node = is.cluster.node)
}
initialize.biclust <-
function (verbose = F, is.cluster.node = F) 
{
    if (!exists("organism") || is.null(organism)) 
        organism <<- params$organism <<- get.organism()
    params$output.dir <<- paste("output/", params$organism, "/", 
        sep = "")
    params$data.dir <<- paste("data/", params$organism, "/", 
        sep = "")
    out.log <<- get.log(prefix = "biclust.log.%03d.txt", output.dir = params$output.dir)
    if (!exists("iter") || is.na(iter)) 
        iter <<- 0
    set.seed(params$rnd.seed)
    params$tmp.prefix <<- paste(params$organism, "_Motif_tmp", 
        sep = "")
    params$state.file <<- paste(params$output.dir, params$organism, 
        sep = "/")
    params$init.state.file <<- paste(params$state.file, ".biclust.init.RData", 
        sep = "")
    params$ps2pdf.bin <<- system("which ps2pdf", intern = T, 
        ignore.stderr = T)
    params$ps2pdf.cmd <<- paste("|", params$ps2pdf.bin, "-")
    attach(params)
    try(detach("global.data"), silent = T)
    load(paste("data/", organism, ".global.data.RData", sep = ""), 
        envir = .GlobalEnv)
    if (exists("global.data")) 
        attach(global.data)
    if (verbose) 
        print.params(params)
    attach(all.iter.params)
    if (verbose) {
        cat.new("ITER-BASED PARAMS:\n")
        print.params(all.iter.params)
    }
    return(params$iter.n.motifs)
}
initializeClusters <-
function () 
{
    return(list())
}
initialize.funcs <-
structure(list(cat.new = function (..., log = get.global("out.log"), 
    file = "", sep = " ", fill = FALSE, labels = NULL, append = FALSE) 
{
    cat(..., file = file, sep = sep, fill = fill, labels = labels, 
        append = append)
    if (file == "" && !is.null(log) && log != "") 
        try(cat(..., file = log, append = T))
}, print.new = function (x, log = get.global("out.log"), ...) 
{
    print(x, ...)
    if (!is.null(log) && log != "") 
        try(capture.output(print(x, ...), file = log, append = T))
}, get.log = function (prefix = "biclust.log.%03d.txt", output.dir = "output/") 
{
    if (exists("out.log")) 
        return(out.log)
    if (!file.exists(output.dir)) 
        try(dir.create(output.dir, recursive = T, showWarnings = F))
    for (i in 1:999) {
        f <- paste(output.dir, sprintf(prefix, as.integer(i)), 
            sep = "")
        if (file.exists(f)) 
            next
        f <- paste(getwd(), f, sep = "/")
        out.log <<- f
        cat.new("Logging all output to", out.log, "\n")
        break
    }
    return(out.log)
}, init = function (is.cluster.node = F) 
{
    initialize.biclust(T, is.cluster.node = is.cluster.node)
}, get.global = function (name) 
{
    return(get(name, .GlobalEnv))
}, get.param = function (name, params = get("params", .GlobalEnv)) 
{
    return(params[[name]])
}, get.iter.based.params = function (name, iter, params = get("params", 
    .GlobalEnv)) 
{
    if (is.null(params[[name]])) 
        name <- paste("iter.", name, sep = "")
    out <- params[[name]][1]
    if (length(params[[name]]) > 1) {
        tmp <- vector()
        tmp[as.integer(names(params[[name]]))] <- params[[name]]
        out <- last.element(tmp[1:iter][!is.na(tmp[1:iter])])
    }
    return(out)
}, get.all.iter.based.params = function (iter = 1, params = get("params", 
    .GlobalEnv), name = NULL, verbose = F) 
{
    if (is.null(name)) {
        out <- list()
        out$iter <- iter
        for (name in names(params)) {
            if (substring(name, 1, 5) != "iter.") 
                next
            val <- get.iter.based.params(name, iter, params)
            real.name <- substring(name, 6)
            out[[real.name]] <- val
        }
        out$r0 <- out$resid.scaling
        out$p0 <- out$motif.scaling
        out$v0 <- out$volume.scaling
        out$q0 <- out$net.scaling
        out$net.q0 <- get.net.q0(params$net.max.weights, iter, 
            out$q0, verbose = verbose)
        return(out)
    }
    else {
        return(get.iter.based.params(name, iter, params))
    }
}, print.params = function (params) 
{
    cat.new("STARTING PARAMETERS:\n\n")
    for (name in sort(names(params))) {
        cat.new(sprintf("%25s", name), "=\t")
        if (class(params[[name]]) == "function") {
            next
        }
        else if (length(params[[name]]) <= 1) {
            cat.new(params[[name]], "\n")
        }
        else {
            if (is.null(names(params[[name]]))) 
                cat.new(which(!is.na(params[[name]])), ": ")
            else cat.new(names(params[[name]])[which(!is.na(params[[name]]))], 
                ": ")
            cat.new(params[[name]][!is.na(params[[name]])], "\n")
        }
    }
}), .Names = c("cat.new", "print.new", "get.log", "init", "get.global", 
"get.param", "get.iter.based.params", "get.all.iter.based.params", 
"print.params"))
last.element <-
function (v) 
{
    return(v[length(v)])
}
load.latest <-
function (kkmax = kmax * 10) 
{
    file <- NULL
    for (i in 1:kkmax) {
        f <- paste(output.dir, "/", organism, ".", sprintf("%03d", 
            as.integer(i)), ".RData", sep = "")
        if (file.exists(f)) 
            file <- f
    }
    if (!is.null(file)) {
        cat.new("Loading saved clusters from", file, "\n")
        load(file)
        cat.new("Cluster stack is now loaded.\n")
        if (is.null(clusterStack$k)) 
            clusterStack$k <- length(clusterStack)
        clusterStack <<- clusterStack
        out.logs <<- out.logs
    }
    else {
        cat.new("Starting run on cluster #1.\n")
    }
}
maxExp <-
function (ratios) 
{
    max.rat <- apply(ratios, 1, max)
    max.rat
}
minExp <-
function (ratios) 
{
    min.rat <- apply(ratios, 1, min)
    min.rat
}
mkBgFile <-
function (bgseqs = NULL, order = 0, bgfname = NULL, input.list = NULL, 
    use.rev.comp = T) 
{
    if (!is.null(input.list)) {
        if (!is.null(bgfname) && !file.exists(bgfname)) {
            tmp <- input.list[2:length(input.list)]
            if (!is.null(bgfname)) 
                write.table(unlist(tmp), row.names = names(tmp), 
                  col.names = paste("#", order, "th order Markov background model"), 
                  quote = F, file = bgfname)
            return(input.list)
        }
    }
    else if (!is.null(bgfname) && file.exists(bgfname)) {
        cat.new("Reading background freqs from", bgfname, "\n")
        tmp <- read.table(bgfname)
        patterns <- as.character(tmp$V1)
        freqs <- as.numeric(tmp$V2)
        order <- max(nchar(patterns)) - 1
        out <- list()
        out$order <- order
        cat.new("BG order was:", order, "\n")
        for (i in 1:length(patterns)) out[[patterns[i]]] <- freqs[i]
        return(out)
    }
    out <- list()
    out$order <- order
    bgseqs <- unique(bgseqs)
    if (use.rev.comp) 
        bgseqs <- unique(c(bgseqs, rev.comp(bgseqs)))
    x <- paste(bgseqs, sep = "", collapse = "")[1]
    for (ord in 0:order) {
        cat.new("Calculating", ord, "th order part of background Markov model from", 
            length(bgseqs), "sequences\n")
        if (use.rev.comp) 
            cat.new("Using reverse-complement too.\n")
        if (ord == 0) {
            bgres <- unlist(strsplit(bgseqs, character(0)), use.names = F)
            bglen <- length(bgres)
            for (res in c("G", "A", "T", "C")) {
                freq <- sum(bgres == res, na.rm = T)/bglen
                out[[res]] <- freq
                cat.new("FREQ:", res, "=", freq, "\n")
            }
            next
        }
        xxx <- character()
        split.size <- 10000
        for (i in 1:(nchar(x)/split.size + 1)) {
            start <- (i - 1) * split.size + 1
            end <- i * split.size + ord
            if (end > nchar(x)) 
                end <- nchar(x)
            xx <- substring(x, start, end)
            len <- nchar(xx)
            xxx <- c(xxx, substring(xx, 1:(len - ord), (1 + ord):len))
            rm(xx)
            if (i%%25 == 0) 
                gc()
        }
        unq <- unique(xxx)
        len <- length(xxx)
        for (i in unq) {
            out[[i]] <- sum(xxx == i)/len
            cat.new("FREQ:", i, "=", out[[i]], "\n")
        }
        rm(xxx)
        gc()
    }
    if (!is.null(bgfname) && !file.exists(bgfname)) {
        cat.new("Writing to file:", bgfname, "\n")
        tmp <- out[2:length(out)]
        write.table(unlist(tmp), row.names = names(tmp), col.names = paste("#", 
            order, "th order Markov background model"), quote = F, 
            file = bgfname)
    }
    invisible(out)
}
mkTempMemeFiles <-
function (sgenes, seqs, fname = "meme.tmp.fst", bgseqs = NULL, 
    bgfname = NULL, bg.list = NULL) 
{
    sgenes <- sgenes[!(is.na(seqs) | is.null(seqs))]
    seqs <- seqs[!(is.na(seqs) | is.null(seqs))]
    lengths <- sum(nchar(seqs)) + length(seqs) * 3
    cat(paste(">", sgenes, "\n", seqs, sep = ""), file = fname, 
        sep = "\n")
    if (!is.null(bgfname) && !file.exists(bgfname)) {
        if (!is.null(bg.list)) 
            mkBgFile(input.list = bg.list, order = bg.list$order, 
                bgfname = bgfname)
        else if (!is.null(bgseqs)) 
            mkBgFile(bgseqs, order = 0, bgfname = bgfname)
    }
    return(lengths)
}
motif.one.cluster <-
function (cluster, force.mcast = F, unlink = T, re.init = F, 
    meme.model = "zoops", verbose = F, n.motifs = NA, mast.return = F) 
{
    if (re.init) 
        initialize.biclust(F)
    if (!exists("upstream.bg") || is.null(upstream.bg) || upstream.bg$order != 
        bg.order) {
        upstream.bg <<- mkBgFile(upstream, order = bg.order, 
            bgfname = paste(data.dir, "/upstream.bg", sep = ""))
    }
    if (force.mcast) 
        use.mcast <- T
    if (is.na(n.motifs)) 
        n.motifs <- all.iter.params$n.motifs
    out <- getBiMotifForCluster(cluster, gene.ids, n.motifs = n.motifs, 
        min.seqs = min.unique.seqs, max.seqs = max.unique.seqs, 
        min.width = min.motif.width, max.width = max.motif.width, 
        do.pal = F, seq.length = seq.length, p.val.cutoff = p.val.cutoff, 
        e.val.cutoff = e.val.cutoff, motif.e.val.cutoff = motif.e.val.cutoff, 
        use.mcast = use.mcast, mast.return = mast.return, meme.model = meme.model, 
        unlink = unlink, verbose = verbose)
    return(out)
}
my.unlink <-
function (files) 
{
    file.remove(files)
}
net.matrix.to.edge.list <-
function (network, edge.names = NULL) 
{
    if (sum(network) == 0) 
        return(matrix(nrow = 1, ncol = 2))
    if (is.null(edge.names)) 
        edge.names <- rownames(network)
    tmp <- which(network > 0, arr.ind = T)
    out <- matrix(nrow = nrow(tmp), ncol = 2)
    out[, 1] <- edge.names[tmp[, "row"]]
    out[, 2] <- edge.names[tmp[, "col"]]
    invisible(out)
}
noChange <-
function (lambdas, lambdaThresh, mode = "names") 
{
    gene.ids <- rownames(lambdas)
    inOrOut <- vector(mode = "numeric", length = length(gene.ids))
    names(inOrOut) <- gene.ids
    noChange <- character()
    for (gene.id in gene.ids) {
        maxLam <- max(lambdas[gene.id, ])
        if (maxLam > lambdaThresh) {
            inOrOut[gene.id] <- 1
        }
        else {
            inOrOut[gene.id] <- 0
            noChange <- c(noChange, gene.id)
        }
    }
    if (mode == "names") 
        return(noChange)
    else return(inOrOut)
}
norm.unit.vector <-
function (xs, ys) 
{
    slope <- (ys[2] - ys[1])/(xs[2] - xs[1])
    norm.slope <- -1/slope
    out <- c(1, norm.slope)
    s <- sqrt(sum(out^2))
    out <- out/s
    return(out)
}
overlap.biclust <-
structure(list(overLapRows = function (rows1, rows2) 
{
    overlap12 <- 0
    rows12 <- sum(rows1 %in% rows2)
    if (rows12 != 0) 
        overlap12 <- rows12/mean(length(rows1), length(rows2))
    return(overlap12)
}, overLapCols = function (cols1, cols2) 
{
    overlap12 <- 0
    cols12 <- sum(cols1 %in% cols2)
    if (cols12 != 0) 
        overlap12 <- cols12/mean(length(cols1), length(cols2))
    return(overlap12)
}, overLapClust = function (clust1, rows) 
{
    return(overLapRows(clust1$rows, rows))
}, cluster.correlation = function (clust1, clust2, ratios = get.global("ratios"), 
    use.cols = c("union", "intersection")) 
{
    if (use.cols == "intersection") 
        cols <- clust1$cols[clust1$cols %in% clust2$cols]
    else if (use.cols == "union") 
        cols <- c(clust1$cols, clust2$cols)
    rat1 <- ratios[clust1$rows, cols]
    prof1 <- apply(rat1, 2, mean, na.rm = T)
    rat2 <- ratios[clust2$rows, cols]
    prof2 <- apply(rat2, 2, mean, na.rm = T)
    cor <- cor(prof1, prof2, use = "pairwise")
    return(cor)
}, compare.cluster.motifs = function (cluster1, cluster2, rev.comp = T, 
    e.v.cutoff = 10) 
{
    m.out1 <- cluster1$motif.out
    m.out2 <- cluster2$motif.out
    if (is.null(m.out1) || is.null(m.out2) || is.null(m.out1$pssms) || 
        is.null(m.out2$pssms) || length(m.out1$pssms) <= 0 || 
        length(m.out2$pssms) <= 0) 
        return(-1)
    cors <- matrix(-1, nrow = length(m.out1$pssms), ncol = length(m.out2$pssms))
    for (i in 1:length(m.out1$pssms)) {
        if (i > ncol(cors)) 
            break
        for (j in i:length(m.out2$pssms)) {
            if (!is.na(m.out1$e.values) && !is.na(m.out2$e.values) && 
                !is.na(m.out1$e.values[i]) && !is.na(m.out2$e.values[j]) && 
                m.out1$e.values[i] <= e.v.cutoff && m.out2$e.values[j] <= 
                e.v.cutoff) {
                cors[i, j] <- compare.pssms(m.out1$pssms[[i]], 
                  m.out2$pssms[[j]], rev.comp = rev.comp)
            }
        }
    }
    cors[cors == 0] <- -1
    return(cors)
}, genes.in.how.many.clusters = function (clusters, genes = get.global("gene.ids"), 
    up.to = clusters$k, use.cols = F) 
{
    out <- integer(length(genes))
    names(out) <- genes
    g <- t(genes)
    if (!use.cols) {
        for (i in 1:up.to) out <- out + (g %in% clusters[[i]]$rows)
    }
    else {
        for (i in 1:up.to) out <- out + (g %in% clusters[[i]]$cols)
    }
    return(out)
}, genes.in.which.clusters = function (clusters, genes = get.global("gene.ids"), 
    min = 1) 
{
    gene.names <- gene.coords$gene.name
    nms <- names(gene.names)
    out <- list()
    if (length(genes) <= 0) 
        return(out)
    for (i in 1:clusters$k) {
        glist <- character()
        for (g in genes) {
            if (!(g %in% nms) && g %in% gene.names) 
                g <- nms[which(gene.names == g)]
            grp <- grep(g, clusters[[i]]$rows, ignore.case = T)
            if (length(grp) > 0) 
                glist <- c(glist, clusters[[i]]$rows[grp])
        }
        if (length(glist) >= min) 
            out[[as.character(i)]] <- glist
    }
    return(out)
}, genes.not.in.any.clusters = function (clusterStack, all.genes = get.global("gene.ids")) 
{
    g <- character()
    for (i in 1:clusterStack$k) g <- c(g, clusterStack[[i]]$rows)
    g <- unique(g)
    out <- all.genes[!all.genes %in% g]
    return(out)
}), .Names = c("overLapRows", "overLapCols", "overLapClust", 
"cluster.correlation", "compare.cluster.motifs", "genes.in.how.many.clusters", 
"genes.in.which.clusters", "genes.not.in.any.clusters"))
overLapClust <-
function (clust1, rows) 
{
    return(overLapRows(clust1$rows, rows))
}
overLapCols <-
function (cols1, cols2) 
{
    overlap12 <- 0
    cols12 <- sum(cols1 %in% cols2)
    if (cols12 != 0) 
        overlap12 <- cols12/mean(length(cols1), length(cols2))
    return(overlap12)
}
overLapRows <-
function (rows1, rows2) 
{
    overlap12 <- 0
    rows12 <- sum(rows1 %in% rows2)
    if (rows12 != 0) 
        overlap12 <- rows12/mean(length(rows1), length(rows2))
    return(overlap12)
}
params <-
structure(list(rnd.seed = 2, n.iters = 100, kmax = 400, stdRatios = TRUE, 
    semirnd.clust.size = 30, ratios.thresh = 0.5, ratios.n.thresh = 10, 
    iter.resid.scaling = 1, r.sig = 0.3, iter.bad.move.temperature = structure(c(0.15, 
    0.148989898989899, 0.147979797979798, 0.146969696969697, 
    0.145959595959596, 0.144949494949495, 0.143939393939394, 
    0.142929292929293, 0.141919191919192, 0.140909090909091, 
    0.13989898989899, 0.138888888888889, 0.137878787878788, 0.136868686868687, 
    0.135858585858586, 0.134848484848485, 0.133838383838384, 
    0.132828282828283, 0.131818181818182, 0.130808080808081, 
    0.12979797979798, 0.128787878787879, 0.127777777777778, 0.126767676767677, 
    0.125757575757576, 0.124747474747475, 0.123737373737374, 
    0.122727272727273, 0.121717171717172, 0.120707070707071, 
    0.11969696969697, 0.118686868686869, 0.117676767676768, 0.116666666666667, 
    0.115656565656566, 0.114646464646465, 0.113636363636364, 
    0.112626262626263, 0.111616161616162, 0.110606060606061, 
    0.10959595959596, 0.108585858585859, 0.107575757575758, 0.106565656565657, 
    0.105555555555556, 0.104545454545455, 0.103535353535354, 
    0.102525252525253, 0.101515151515152, 0.10050505050505, 0.0994949494949495, 
    0.0984848484848485, 0.0974747474747475, 0.0964646464646465, 
    0.0954545454545455, 0.0944444444444444, 0.0934343434343434, 
    0.0924242424242424, 0.0914141414141414, 0.0904040404040404, 
    0.0893939393939394, 0.0883838383838384, 0.0873737373737374, 
    0.0863636363636364, 0.0853535353535353, 0.0843434343434343, 
    0.0833333333333333, 0.0823232323232323, 0.0813131313131313, 
    0.0803030303030303, 0.0792929292929293, 0.0782828282828283, 
    0.0772727272727273, 0.0762626262626263, 0.0752525252525253, 
    0.0742424242424242, 0.0732323232323232, 0.0722222222222222, 
    0.0712121212121212, 0.0702020202020202, 0.0691919191919192, 
    0.0681818181818182, 0.0671717171717172, 0.0661616161616162, 
    0.0651515151515151, 0.0641414141414141, 0.0631313131313131, 
    0.0621212121212121, 0.0611111111111111, 0.0601010101010101, 
    0.0590909090909091, 0.0580808080808081, 0.0570707070707071, 
    0.0560606060606061, 0.055050505050505, 0.054040404040404, 
    0.053030303030303, 0.052020202020202, 0.051010101010101, 
    0.05), .Names = c("1", "2", "3", "4", "5", "6", "7", "8", 
    "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", 
    "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", 
    "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", 
    "39", "40", "41", "42", "43", "44", "45", "46", "47", "48", 
    "49", "50", "51", "52", "53", "54", "55", "56", "57", "58", 
    "59", "60", "61", "62", "63", "64", "65", "66", "67", "68", 
    "69", "70", "71", "72", "73", "74", "75", "76", "77", "78", 
    "79", "80", "81", "82", "83", "84", "85", "86", "87", "88", 
    "89", "90", "91", "92", "93", "94", "95", "96", "97", "98", 
    "99", "100")), cluster.row.floor = 3, mean.clusters.per.gene = 2, 
    expected.cluster.rows = 30, weight.fixed = 0.5, max.add.rows.per.iter = 3, 
    max.remove.rows.per.iter = 3, max.col.moves.per.cluster = 5, 
    allow.anticor = FALSE, iter.motif.scaling = structure(c(0, 
    0, 0.0258620689655172, 0.0517241379310345, 0.0775862068965517, 
    0.103448275862069, 0.129310344827586, 0.155172413793103, 
    0.181034482758621, 0.206896551724138, 0.232758620689655, 
    0.258620689655172, 0.28448275862069, 0.310344827586207, 0.336206896551724, 
    0.362068965517241, 0.387931034482759, 0.413793103448276, 
    0.439655172413793, 0.46551724137931, 0.491379310344828, 0.517241379310345, 
    0.543103448275862, 0.568965517241379, 0.594827586206897, 
    0.620689655172414, 0.646551724137931, 0.672413793103448, 
    0.698275862068966, 0.724137931034483, 0.75, 0.75, 0.739795918367347, 
    0.729591836734694, 0.719387755102041, 0.709183673469388, 
    0.698979591836735, 0.688775510204082, 0.678571428571429, 
    0.668367346938776, 0.658163265306122, 0.647959183673469, 
    0.637755102040816, 0.627551020408163, 0.61734693877551, 0.607142857142857, 
    0.596938775510204, 0.586734693877551, 0.576530612244898, 
    0.566326530612245, 0.556122448979592, 0.545918367346939, 
    0.535714285714286, 0.525510204081633, 0.51530612244898, 0.505102040816326, 
    0.494897959183673, 0.48469387755102, 0.474489795918367, 0.464285714285714, 
    0.454081632653061, 0.443877551020408, 0.433673469387755, 
    0.423469387755102, 0.413265306122449, 0.403061224489796, 
    0.392857142857143, 0.38265306122449, 0.372448979591837, 0.362244897959184, 
    0.352040816326531, 0.341836734693878, 0.331632653061225, 
    0.321428571428571, 0.311224489795918, 0.301020408163265, 
    0.290816326530612, 0.280612244897959, 0.270408163265306, 
    0.260204081632653, 0.25), .Names = c("1", "5", "6", "7", 
    "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", 
    "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", 
    "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", 
    "38", "39", "40", "41", "42", "43", "44", "45", "46", "47", 
    "48", "49", "50", "51", "52", "53", "54", "55", "56", "57", 
    "58", "59", "60", "61", "62", "63", "64", "65", "66", "67", 
    "68", "69", "70", "71", "72", "73", "74", "75", "76", "77", 
    "78", "79", "80", "81", "82", "83", "84")), iter.min.unique.seqs = 3, 
    iter.max.unique.seqs = 200, iter.n.motifs = structure(c(1, 
    2, 2), .Names = c("1", "25", "40")), iter.min.motif.width = 9, 
    iter.max.motif.width = 25, iter.seq.length = 250, iter.bg.order = 3, 
    iter.e.val.cutoff = structure(c(999999, 99999, 9999), .Names = c("1", 
    "25", "40")), iter.p.val.cutoff = structure(c(0.9, 0.9), .Names = c("1", 
    "65")), iter.motif.e.val.cutoff = structure(c(999999, 9999
    ), .Names = c("1", "40")), iter.meme.model = structure("zoops", .Names = "1"), 
    net.max.weights = structure(c(2, 0.8, 0, 0.1, 0, 0, 0, 0, 
    0, 0, 0.1), .Names = c("operons", "met", "prolinks.RS", "prolinks.PP", 
    "prolinks.GN", "prolinks.GC", "COG.code", "cond.sims", "predictome.CP", 
    "predictome.GF", "predictome.PP")), iter.net.scaling = structure(c(0.1, 
    0.125, 0.15, 0.175, 0.2, 0.2, 0.196326530612245, 0.19265306122449, 
    0.188979591836735, 0.18530612244898, 0.181632653061225, 0.177959183673469, 
    0.174285714285714, 0.170612244897959, 0.166938775510204, 
    0.163265306122449, 0.159591836734694, 0.155918367346939, 
    0.152244897959184, 0.148571428571429, 0.144897959183673, 
    0.141224489795918, 0.137551020408163, 0.133877551020408, 
    0.130204081632653, 0.126530612244898, 0.122857142857143, 
    0.119183673469388, 0.115510204081633, 0.111836734693878, 
    0.108163265306122, 0.104489795918367, 0.100816326530612, 
    0.0971428571428571, 0.0934693877551021, 0.0897959183673469, 
    0.0861224489795918, 0.0824489795918367, 0.0787755102040816, 
    0.0751020408163265, 0.0714285714285714, 0.0677551020408163, 
    0.0640816326530612, 0.0604081632653061, 0.056734693877551, 
    0.0530612244897959, 0.0493877551020408, 0.0457142857142857, 
    0.0420408163265306, 0.0383673469387755, 0.0346938775510204, 
    0.0310204081632653, 0.0273469387755102, 0.0236734693877551, 
    0.02), .Names = c("1", "2", "3", "4", "5", "11", "12", "13", 
    "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", 
    "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", 
    "34", "35", "36", "37", "38", "39", "40", "41", "42", "43", 
    "44", "45", "46", "47", "48", "49", "50", "51", "52", "53", 
    "54", "55", "56", "57", "58", "59", "60")), operon.nets = c("prolinks.GC", 
    "predictome.CP"), pdf.n.best = 300, meme.cmd = "./progs/meme", 
    mast.cmd = "./progs/mast", tmp.dir = "/tmp/RtmpL1hRjk", data.dir = "data/halo/", 
    ps2pdf.bin = "/usr/bin/ps2pdf", ps2pdf.cmd = "| /usr/bin/ps2pdf - ", 
    operon.net = "operons", plot.it = FALSE), .Names = c("rnd.seed", 
"n.iters", "kmax", "stdRatios", "semirnd.clust.size", "ratios.thresh", 
"ratios.n.thresh", "iter.resid.scaling", "r.sig", "iter.bad.move.temperature", 
"cluster.row.floor", "mean.clusters.per.gene", "expected.cluster.rows", 
"weight.fixed", "max.add.rows.per.iter", "max.remove.rows.per.iter", 
"max.col.moves.per.cluster", "allow.anticor", "iter.motif.scaling", 
"iter.min.unique.seqs", "iter.max.unique.seqs", "iter.n.motifs", 
"iter.min.motif.width", "iter.max.motif.width", "iter.seq.length", 
"iter.bg.order", "iter.e.val.cutoff", "iter.p.val.cutoff", "iter.motif.e.val.cutoff", 
"iter.meme.model", "net.max.weights", "iter.net.scaling", "operon.nets", 
"pdf.n.best", "meme.cmd", "mast.cmd", "tmp.dir", "data.dir", 
"ps2pdf.bin", "ps2pdf.cmd", "operon.net", "plot.it"))
parseMastDiagram <-
function (diagram, motif.widths) 
{
    motif <- posn <- width <- numeric()
    arr <- strsplit(diagram, "_")[[1]]
    motif.inds <- grep("[\\]\\>]", arr, perl = T)
    arr <- as.vector(sapply(arr, function(i) sub("\\[", "", i)))
    arr <- as.vector(sapply(arr, function(i) sub("<", "", i)))
    arr <- as.vector(sapply(arr, function(i) sub("\\]", "", i)))
    arr <- as.numeric(as.vector(sapply(arr, function(i) sub(">", 
        "", i))))
    loc <- 1
    ind <- i <- 1
    while (i <= length(arr)) {
        if (!(i %in% motif.inds)) 
            loc <- loc + arr[i]
        else {
            motif[ind] <- arr[i]
            posn[ind] <- loc
            width[ind] <- motif.widths[abs(arr[i])]
            loc <- loc + width[ind]
            ind <- ind + 1
        }
        i <- i + 1
    }
    tmp <- cbind(motif, posn, width)
    if (nrow(tmp) <= 0) 
        return(NULL)
    return(tmp)
}
plotCluster <-
function (cluster, ratios = get.global("ratios"), iter = get.global("iter"), 
    postscript.file = NULL, regulators = NULL, cond.labels = F) 
{
    if (!is.null(postscript.file)) 
        postscript(postscript.file)
    k <- cluster$k
    titl <- paste("Cluster:", k, "; resid:", sprintf("%.3f", 
        as.numeric(cluster$resid)), "; genes:", length(cluster$rows), 
        "; conds:", length(cluster$cols), "; iter:", iter)
    range.r <- range(ratios[cluster$rows, cluster$cols], na.rm = T)
    if (cluster$ncols < 100) 
        range.r[1] <- range.r[1] * 1.5
    cols.b <- colnames(ratios)[colnames(ratios) %in% cluster$cols]
    old.pars <- par()
    if (!cond.labels) 
        par(mar = rep(2, 4), mgp = c(3, 1, 0) * 0.5)
    plot(1:length(cols.b), ratios[cluster$rows[1], cols.b], ylim = range.r, 
        type = "l", main = titl, cex.main = 0.7, xlab = NA, ylab = NA, 
        cex.lab = 0.7, cex.sub = 0.7, cex.axis = 0.7, lwd = 1)
    colmap <- rainbow(cluster$nrows)
    for (i in 1:cluster$nrows) {
        lines(1:length(cols.b), ratios[cluster$rows[i], cols.b], 
            col = colmap[i], lwd = 1)
    }
    if (cond.labels) {
        tmp.y <- (1:cluster$ncols/1:cluster$ncols) * range.r[1] * 
            0.85
        text(1:cluster$ncols, tmp.y, cols.b, srt = 90, cex = 0.2)
    }
    if (!is.null(regulators)) {
        legend <- character()
        cols <- ltys <- lwds <- integer()
        col.i <- 1
        for (i in 1:length(regulators)) {
            if (length(grep(".m[ai][nx]", regulators[i])) > 0) {
                splt <- strsplit(regulators[i], "\\.")[[1]]
                lines(1:length(cols.b), apply(ratios[splt[1:2], 
                  cols.b], 2, splt[3]), col = col.i, lty = 1, 
                  lwd = 3)
                legend <- c(legend, regulators[i])
                cols <- c(cols, col.i)
                ltys <- c(ltys, 1)
                lwds <- c(lwds, 3)
                col.i <- col.i + 1
            }
            else {
                lines(1:length(cols.b), ratios[regulators[i], 
                  cols.b], col = col.i, lty = 1, lwd = 3)
                legend <- c(legend, regulators[i])
                cols <- c(cols, col.i)
                ltys <- c(ltys, 1)
                lwds <- c(lwds, 3)
                col.i <- col.i + 1
            }
        }
        legend(length(cols.b) * 2/3, range.r[2], legend, lty = ltys, 
            lwd = lwds, col = cols, cex = 0.7)
    }
    par(old.pars)
    if (!is.null(postscript.file)) 
        graphics.off()
}
plotCluster.all.conds <-
function (cluster, ratios = get.global("ratios"), iter = get.global("iter"), 
    postscript.file = NULL, plot.resids = T, regulators = NULL, 
    cond.labels = F) 
{
    if (!is.null(postscript.file)) 
        postscript(postscript.file, paper = "letter")
    k <- cluster$k
    titl <- paste("Cluster:", k, "; resid:", sprintf("%.3f", 
        as.numeric(cluster$resid)), "; genes:", length(cluster$rows), 
        "; conds:", length(cluster$cols), "; iter:", iter)
    range.r <- range(ratios[cluster$rows, ], na.rm = T)
    cols.b <- colnames(ratios)[colnames(ratios) %in% cluster$cols]
    cols.b <- c(cols.b, colnames(ratios)[!colnames(ratios) %in% 
        cluster$cols])
    len.b <- length(cols.b)
    if (len.b < 100) 
        range.r[1] <- range.r[1] * 1.5
    old.pars <- par()
    if (!cond.labels) 
        par(mar = rep(2, 4), mgp = c(3, 1, 0) * 0.5)
    rats <- ratios[cluster$rows, cols.b]
    plot(1:len.b, rats[cluster$rows[1], cols.b], ylim = range.r, 
        type = "l", col = 1, main = titl, cex.main = 0.7, xlab = NA, 
        ylab = NA, cex.lab = 0.7, cex.sub = 0.7, cex.axis = 0.7)
    if (plot.resids) {
        resids <- residOneClust.by.col(rats)
        resids <- resids/max(resids) * max(rats, na.rm = T)
        poly.x <- c(1:length(cols.b), length(cols.b):1)
        poly.y <- c(resids, rev(-resids))
        polygon(poly.x, poly.y, border = NA, col = "lightgray")
    }
    colmap <- rainbow(cluster$nrows)
    for (i in 1:cluster$nrows) {
        lines(1:len.b, rats[cluster$rows[i], cols.b], col = colmap[i], 
            type = "l", lwd = 1)
    }
    cut.x <- c(cluster$ncols + 0.5, cluster$ncols + 0.5)
    lines(cut.x, range.r, col = 2, lwd = 3, lty = 2)
    if (!is.null(regulators)) {
        legend <- character()
        cols <- ltys <- lwds <- integer()
        col.i <- 1
        for (i in 1:length(regulators)) {
            if (length(grep(".m[ai][nx]", regulators[i])) > 0) {
                splt <- strsplit(regulators[i], "\\.")[[1]]
                lines(1:length(cols.b), apply(ratios[splt[1:2], 
                  cols.b], 2, splt[3]), col = col.i, lty = 1, 
                  lwd = 3)
                legend <- c(legend, regulators[i])
                cols <- c(cols, col.i)
                ltys <- c(ltys, 1)
                lwds <- c(lwds, 3)
                col.i <- col.i + 1
            }
            else {
                lines(1:length(cols.b), ratios[regulators[i], 
                  cols.b], col = col.i, lty = 1, lwd = 3)
                legend <- c(legend, regulators[i])
                cols <- c(cols, col.i)
                ltys <- c(ltys, 1)
                lwds <- c(lwds, 3)
                col.i <- col.i + 1
            }
        }
        legend(length(cols.b) * 2/3, range.r[2], legend, lty = ltys, 
            lwd = lwds, col = cols, cex = 0.7)
    }
    if (cond.labels) {
        tmp.y <- (1:len.b/1:len.b) * range.r[1] * 0.85
        text(1:len.b, tmp.y, cols.b, srt = 90, cex = 0.2)
    }
    par(old.pars)
    if (!is.null(postscript.file)) 
        graphics.off()
}
plotCluster.motif <-
function (cluster, iter = get.global("iter"), ratios = get.global("ratios"), 
    seqs = get.global("opUpstream"), postscript.file = NULL, 
    min.e.val = 9998, both.motifs = F, all.conds = T, networks = get.global("networks"), 
    short.names = T, shade = T, print.pvalues = F, regulators = NULL, 
    separate.pages = F, force.labels = F) 
{
    if (!is.null(postscript.file)) 
        postscript(postscript.file, paper = "letter", horiz = F)
    if (!separate.pages) {
        split.screen(c(2, 2))
        screen(screen <- 1)
    }
    k <- cluster$k
    if (!is.null(ratios)) {
        if (all.conds || !both.motifs) {
            plotCluster.all.conds(cluster, ratios, postscript = NULL, 
                regulators = regulators, cond.labels = separate.pages | 
                  force.labels, iter = iter)
            if (!separate.pages) 
                screen(screen <- screen + 1)
        }
        if (!all.conds || !both.motifs) {
            plotCluster(cluster, ratios, postscript = NULL, regulators = regulators, 
                cond.labels = separate.pages | force.labels, 
                iter = iter)
            if (!separate.pages) 
                screen(screen <- screen + 1)
        }
    }
    if (!is.null(regulators)) {
        cluster$rows <- unique(c(cluster$rows, regulators))
        cluster$nrows <- length(cluster$rows)
    }
    if (!is.null(cluster$e.val) && !is.null(cluster$motif.out)) {
        if (min(cluster$e.val, na.rm = T) < min.e.val && !is.na(cluster$p.clust)) {
            if (both.motifs) {
                meme.out <- cluster$motif.out$meme.out$nonpal
                if (!is.null(meme.out)) {
                  memeOut.info <- getMemeMotifInfo(meme.out)[[1]]
                  viewPssm(memeOut.info$pssm, e.val = memeOut.info$e.value, 
                    mode = "scale", separate.pages = separate.pages)
                  if (!separate.pages) 
                    screen(screen <- screen + 1)
                }
                meme.out <- cluster$motif.out$meme.out$pal
                if (!is.null(meme.out)) {
                  memeOut.info <- getMemeMotifInfo(meme.out)[[1]]
                  viewPssm(memeOut.info$pssm, e.val = memeOut.info$e.value, 
                    p.clust = cluster$p.clust, mode = "scale", 
                    separate.pages = separate.pages)
                  if (!separate.pages) 
                    screen(screen <- screen + 1)
                }
            }
            else {
                if (!is.null(cluster$motif.out$pssms)) {
                  pssm <- cluster$motif.out$pssms
                  if (length(pssm) > 0) {
                    if (!separate.pages) {
                      subscreens <- split.screen(c(2, 2), screen = screen)
                      subscreen <- 1
                    }
                    else {
                      par(mfrow = c(2, 2))
                    }
                    for (ppp in 1:length(pssm)) {
                      viewPssm(pssm[[ppp]], e.val = cluster$motif.out$e.values[ppp], 
                        mode = "scale", mot.ind = ppp, separate.pages = separate.pages)
                      if (!separate.pages) 
                        subscreen <- subscreen + 1
                      if (!separate.pages) 
                        screen(subscreens[subscreen])
                    }
                    if (separate.pages) 
                      par(mfrow = c(1, 1))
                  }
                }
            }
        }
        tmp <- plot.cluster.network(cluster, networks, gene.ids)
        if (!separate.pages) 
            screen(screen <- screen + 1)
        if (!is.null(seqs)) 
            plotClusterMotifPositions(cluster, seqs, short = short.names, 
                shade = shade, colors = tmp$vizmap$ncol)
    }
    else if (!is.null(networks)) {
        tmp <- plot.cluster.network(cluster, networks, gene.ids)
        if (!separate.pages) 
            screen(screen <- screen + 1)
        if (!is.null(seqs)) 
            plotClusterMotifPositions(cluster, seqs, short = short.names, 
                shade = shade)
    }
    if (print.pvalues && !is.null(networks)) 
        print.net.pvalues(k, clusterStack, networks, gene.ids)
    if (!separate.pages) 
        close.screen(all = TRUE)
    if (!is.null(postscript.file)) 
        graphics.off()
}
plotClusterMotifPositions <-
function (cluster, seqs = get("opUpstream", .GlobalEnv), short.names = T, 
    shade = T, p.val.shade.cutoff = NA, colors = NULL) 
{
    rows <- cluster$rows
    motif.out <- cluster$motif.out
    is.dup.seq <- get.cluster.dup.seqs(cluster)
    p.clust = cluster$p.clust
    no.motif <- FALSE
    p.values <- motif.widths <- pssm <- diagrams <- NULL
    if (!is.null(motif.out) && !is.null(motif.out$diagrams)) {
        p.values <- motif.out$p.values[rows]
        motif.widths <- sapply(motif.out$pssms, nrow, simplify = T)
        pssm <- motif.out$pssms
        diagrams <- sapply(rows, function(i) motif.out$diagrams[i])
    }
    else {
        no.motif <- TRUE
        diagrams <- character(length(rows))
        p.values <- numeric(length(rows))
        motif.widths <- 0
    }
    seqs <- seqs[rows]
    names(seqs) <- rows
    names(diagrams) <- rows
    names(p.values) <- rows
    if (!is.null(diagrams) && length(diagrams) != length(seqs)) 
        cat.new("WARNING -- diagrams and seqs not the same length!")
    seq.lengths <- sapply(seqs, function(i) nchar(i))
    if (any(seq.lengths > median(seq.lengths))) {
        seqs <- subSeq(seqs, len = median(seq.lengths), direction = "back")
        seq.lengths <- sapply(seqs, function(i) nchar(i))
    }
    if (is.na(p.val.shade.cutoff)) 
        p.val.shade.cutoff <- p.val.cutoff
    if (use.mcast) 
        p.val.shade.cutoff <- 0.001
    maxlen <- max(seq.lengths)
    g.names <- names(diagrams)
    n.genes <- length(diagrams)
    inds <- integer()
    sort.by <- ""
    if (no.motif) 
        sort.by <- "gene.name"
    if (sort.by == "gene.name") 
        inds <- sort(g.names, decreas = T, index = T)$ix
    else inds <- sort(p.values[g.names], decreas = T, index = T)$ix
    x.range <- c(-maxlen * 0.15, maxlen * 1.08)
    y.range <- c(0.5, length(diagrams) + 1)
    old.pars <- par()
    par(mar = rep(0.5, 4), mgp = c(3, 1, 0) * 0.5)
    plot(x.range, y.range, xlab = "sequence position", ylab = "sequence", 
        type = "n", axes = F, tck = 0.01, cex.lab = 0.2, cex.sub = 0.2, 
        cex.axis = 0.2)
    axis(side = 1, pos = 0.6, cex.lab = 0.5, cex.sub = 0.5, cex.axis = 0.5, 
        tck = 0.01, mgp = c(0.1, 0.1, 0.1))
    mots.used <- numeric()
    colmap <- rainbow(length(diagrams))
    if (is.list(motif.widths)) {
        if (length(motif.widths) <= 0) 
            motif.widths <- 0
        else {
            for (i in 1:length(motif.widths)) if (is.null(motif.widths[[i]])) 
                motif.widths[[i]] <- 0
            motif.widths <- unlist(motif.widths)
        }
    }
    for (j in 1:length(diagrams)) {
        jj <- inds[j]
        cur.gene <- names(diagrams)[jj]
        if (!is.na(diagrams[jj]) && !is.null(diagrams[jj]) && 
            diagrams[jj] != "") {
            parsed <- parseMastDiagram(diagrams[jj], motif.widths)
            if (!is.null(parsed) && !is.na(parsed)) {
                mot.info <- motif.out$mast.info[[cur.gene]]
                mots <- parsed[, "motif"]
                starts <- parsed[, "posn"]
                widths <- parsed[, "width"]
                for (i in 1:length(mots)) {
                  mot <- mots[i]
                  if (is.na(mot)) 
                    next
                  mots.used <- append(mots.used, abs(mot))
                  start <- starts[i]
                  end <- start + widths[i]
                  col <- abs(mot) + 1
                  if (shade) {
                    if (!is.null(mot.info)) 
                      p.val <- mot.info$pvals[i]
                    else p.val <- 1e-05
                    if (is.na(p.val) || p.val > p.val.shade.cutoff || 
                      p.val > 1) 
                      next
                    else if (p.val <= 0) 
                      p.val <- 1e-05
                    p.val <- log10(p.val)
                    p.min <- -5
                    p.max <- log10(p.val.shade.cutoff)
                    col <- col2rgb(palette()[col])/255
                    col[col > 0] <- 1
                    if (p.val < 10) 
                      col[col == 0] <- (p.val - p.min)/(p.max - 
                        p.min)
                    else col[col == 0] <- 0.99
                    col[col < 0] <- 0
                    col[col > 1] <- 1
                    col <- rgb(col["red", 1], col["green", 1], 
                      col["blue", 1])
                  }
                  if (!is.null(mot.info)) {
                    if (mot > 0) 
                      rect(start, j + 0.01, end, j + 0.3, col = col, 
                        border = col)
                    else if (mot < 0) 
                      rect(start, j - 0.3, end, j - 0.01, col = col, 
                        border = col)
                  }
                  else {
                    if (mot > 0) 
                      rect(start, j + 0.01, end, j + 0.3, border = col)
                    else if (mot < 0) 
                      rect(start, j - 0.3, end, j - 0.01, border = col)
                  }
                }
            }
        }
        if (!is.null(g.names)) {
            label <- g.names[jj]
            if (!is.null(colors)) 
                rect(-maxlen * 0.195, j - 0.18, -5, j + 0.18, 
                  col = NA, border = colors[label], lwd = 3)
        }
        lwd <- 3
        if (n.genes > 40) 
            lwd <- 1
        else if (n.genes > 20) 
            lwd <- 2
        if (jj == 1) 
            lines(c(1, seq.lengths[jj]), c(j, j), lwd = lwd)
        else lines(c(1, seq.lengths[jj]), c(j, j), lwd = lwd, 
            col = colmap[jj])
        if (!is.null(g.names)) {
            label <- g.names[jj]
            if (exists("all.tfs") && label %in% all.tfs) 
                text(-maxlen * 0.2, j, labels = label, cex = 0.7, 
                  adj = c(0, 0.5), col = "tomato3")
            else text(-maxlen * 0.2, j, labels = label, cex = 0.7, 
                adj = c(0, 0.5))
            if (short.names) {
                g.name <- g.names[jj]
                if (exists("gene.coords") && !is.null(gene.coords$gene.func)) 
                  g.name <- gene.coords$gene.func[g.name]
                if (exists("gene.coords") && !is.null(gene.coords$gene.name) && 
                  gene.coords$gene.name[g.names[jj]] != "-" && 
                  !is.na(gene.coords$gene.name[g.names[jj]]) && 
                  toupper(gene.coords$gene.name[g.names[jj]]) != 
                    toupper(g.name)) 
                  g.name <- paste(gene.coords$gene.name[g.names[jj]], 
                    ": ", g.name, sep = "")
                if (!is.na(g.name) && toupper(g.name) != toupper(g.names[jj])) 
                  text(-maxlen * 0.2, j + 0.3, labels = g.name, 
                    cex = 0.7, adj = c(0, 0.5))
            }
            text(maxlen * 1.07, j, labels = prettyNum(p.values[g.names[jj]], 
                digits = 2), cex = 0.7, xpd = NA)
        }
        else text(-maxlen * 0.1, jj, labels = as.character(j), 
            cex = 0.7)
    }
    text(maxlen * 1.05, length(diagrams) + 0.9, labels = "log10(P)", 
        cex = 0.7)
    mots.used <- sort(unique(mots.used))
    if (length(mots.used) > 1) {
        text(-maxlen * 0.2, length(diagrams) + 0.9, "Motif legend:", 
            xpd = NA, adj = c(0, 0.5))
        for (j in 1:length(mots.used)) text((j + 1) * maxlen * 
            0.03, length(diagrams) + 0.9, as.character(mots.used[j]), 
            col = mots.used[j] + 1, xpd = NA, adj = c(0, 0.5))
    }
    if (!is.na(p.clust)) 
        text(maxlen * 0.3, length(diagrams) + 0.9, paste("log10(P.clust)=", 
            sprintf("%.3f", as.numeric(p.clust))), xpd = NA, 
            cex = 0.7)
    n.unique.seqs <- sum(!is.dup.seq)
    text(maxlen * 0.7, length(diagrams) + 0.9, paste(length(seqs), 
        "sequences;", n.unique.seqs, "unique"), xpd = NA, cex = 0.7)
    par(old.pars)
}
plot.cluster.network <-
function (cluster, networks = get("networks", .GlobalEnv), gene.ids = get("gene.ids", 
    .GlobalEnv), common.names = NULL, edge.pvals = F, animate = 0, 
    layout = c("circle", "spring"), skip = "COG.code") 
{
    if (is.null(cluster$net.p.old)) 
        cluster$net.p.old <- subnet.p.values(cluster$rows, networks, 
            gene.ids, edge = T)
    l.tmp <- cluster$net.p.old
    genes <- cluster$rows
    names(genes) <- NULL
    edges <- NULL
    nodes <- character()
    edge.types <- character()
    type.lookup <- 1:length(names(l.tmp))
    names(type.lookup) <- names(l.tmp)
    nodes <- genes
    gene.inds <- which(gene.ids %in% genes)
    for (ii in names(l.tmp)) {
        if (ii == "cond.sims" || ii %in% skip) 
            next
        network <- networks[[ii]]
        tmp <- network[, 1] %in% genes
        if (sum(tmp) <= 0) 
            next
        tmp <- network[tmp, ]
        if (is.vector(tmp)) 
            tmp <- matrix(tmp, ncol = 2)
        tmp2 <- tmp[, 2] %in% genes
        if (sum(tmp2) <= 0) 
            next
        tmp <- tmp[tmp2, ]
        edges2 <- uniquify.edge.list(tmp)
        nodes <- c(nodes, genes)
        match <- sapply(gene.ids[gene.inds], function(i) which(genes == 
            i))
        edges2 <- matrix(match[edges2], ncol = 2)
        edges <- rbind(edges, edges2)
        edge.types <- c(edge.types, rep(ii, nrow(edges2)))
    }
    if (is.null(edges) || nrow(edges) <= 1) {
        edges <- NULL
    }
    else {
        edges <- edges[2:nrow(edges), ]
        if (!is.matrix(edges)) 
            edges <- matrix(edges, ncol = 2)
    }
    nodes <- unique(nodes)
    if (!is.null(common.names)) {
        new.nodes <- common.names[nodes]
        bad.nodes <- which(is.na(new.nodes) | new.nodes == "-")
        new.nodes[bad.nodes] <- nodes[bad.nodes]
        nodes <- new.nodes
    }
    ecols <- rainbow(length(unique(edge.types)))
    names(ecols) <- unique(edge.types)
    vizmap <- get.vizmap(nodes, edges, nsize = 2, ncol = "darkgrey", 
        nfill = "white", lsize = 0.7, lcol = "black", ecol = ecols[edge.types])
    tmp.lett <- 1:26
    names(tmp.lett) <- LETTERS
    coo <- gene.coords$gene.code[nodes]
    tmp <- unique(tmp.lett[coo])
    names(tmp) <- names(tmp.lett[coo][!duplicated(tmp.lett[coo])])
    cols <- rainbow(length(tmp))
    names(cols) <- names(tmp)
    cols <- cols[names(tmp.lett[coo])]
    cols[is.na(names(cols))] <- "darkgrey"
    names(cols) <- genes
    vizmap$ncol <- vizmap$nfill <- cols
    vizmap$lcol <- rep("black", length(nodes))
    if (exists("all.tfs")) 
        vizmap$lcol[which(nodes %in% all.tfs)] <- "tomato3"
    if (layout == "circle") 
        node.coords <- CircleGraphLayout(nodes, rad = 5)
    else if (layout == "spring") 
        node.coords <- SpringEmbeddedLayoutAlgorithm(nodes, edges, 
            node.coords = NULL, vizmap = vizmap, plot = animate, 
            stop.disp = 0.25, k.fac = 0.1)
    old.pars <- par()
    par(mar = rep(0.5, 4), mgp = c(3, 1, 0) * 0.5)
    plot.network(nodes, edges, node.coords, vizmap = vizmap, 
        cex.main = 0.01, axes = F, xlab = NA, ylab = NA, tck = 0.01, 
        cex.lab = 0.2, cex.sub = 0.2, cex.axis = 0.2)
    if (!is.null(edges)) {
        e.types <- unique(names(vizmap$ecol))
        tmp <- legend(0, 2, e.types, xjust = 0.5, lty = 1, cex = 0.7, 
            col = vizmap$ecol[e.types])
        tmp.y <- tmp$rect$top - tmp$rect$h - 0.1
    }
    else {
        tmp.y <- 0
    }
    node.names <- names(cols)
    col <- cols[node.names]
    names(col) <- gene.coords$gene.code[node.names]
    tmp <- unique(gene.coords$gene.code[node.names])
    tmp <- tmp[tmp != "-"]
    if (length(tmp) > 0) 
        legend(0, tmp.y, tmp, xjust = 0.5, lty = 0, pch = 19, 
            cex = 0.7, col = col[tmp])
    invisible(list(nodes = nodes, edges = edges, types = edge.types, 
        coords = node.coords, vizmap = vizmap))
}
plot.clusters <-
function (clusters) 
{
    p0 <<- 1
    plotClusterStack.motif(clusters, ratios, opUpstream, sort = "e.norm", 
        shade = T, n.best = pdf.n.best, networks = networks, 
        postscript = paste(ps2pdf.cmd, " ", output.dir, "/zzz_enorm.pdf", 
            sep = ""))
}
plotClusterStack <-
function (clust = get("clusterStack", .GlobalEnv), ratios = get("ratios", 
    .GlobalEnv), cNum = "all", postscript.file = NULL, sort.em = T) 
{
    if (!is.null(postscript.file)) 
        postscript(postscript.file, paper = "letter")
    if (cNum == "all") {
        sorted <- list()
        sorted$ix <- 1:clust$k
        if (sort.em) {
            resids <- unlist(lapply(clust, function(i) i$resid), 
                use.names = F)
            sorted <- sort(resids, index.return = T)
        }
        for (kk in 1:clust$k) {
            k <- sorted$ix[kk]
            plotCluster(clust[[k]], ratios, postscript = NULL)
            if (is.null(postscript.file)) {
                tmp.menu <- menu(c(" next plot", " quit"))
                if (tmp.menu == 0 || tmp.menu == 2) 
                  break
            }
        }
    }
    else {
        k <- cNum
        plotCluster(clust[[k]], ratios, postscript = NULL)
    }
    if (!is.null(postscript.file)) 
        graphics.off()
}
plotClusterStack.all.conds <-
function (clust = get.global("clusterStack"), ratios = get.global("ratios"), 
    cNum = "all", postscript.file = NULL, sort.em = T, no.menu = F) 
{
    if (!is.null(postscript.file)) 
        postscript(postscript.file, paper = "letter")
    if (cNum == "all") {
        sorted <- list()
        sorted$ix <- 1:clust$k
        if (sort.em) {
            resids <- unlist(lapply(clust, function(i) i$resid), 
                use.names = F)
            sorted <- sort(resids, index.return = T)
        }
        for (kk in 1:clust$k) {
            k <- sorted$ix[kk]
            plotCluster.all.conds(clust[[k]], ratios, postscript = NULL)
            cat.new(clust[[k]]$rows, "\n", gene.coords$gene.name[clust[[k]]$rows], 
                "\n")
            if (is.null(postscript.file) || !no.menu) {
                tmp.menu <- menu(c(" next plot", " quit"))
                if (tmp.menu == 0 || tmp.menu == 2) 
                  break
            }
        }
    }
    else {
        k <- cNum
        plotCluster.all.conds(clust[[k]], ratios, postscript = NULL)
    }
    if (!is.null(postscript.file)) 
        graphics.off()
}
plotClusterStack.motif <-
function (clusterStack = get.global("clusterStack"), ratios = get.global("ratios"), 
    seqs = get.global("opUpstream"), postscript.file = NULL, 
    iter = get.global("iter"), sort.by = c("e.norm", "p.clust", 
        "e.val", "resid"), min.e.val = 9998, ks = 1:clusterStack$k, 
    short.names = T, both.motifs = F, networks = get.global("networks"), 
    shade = T, n.best = 100) 
{
    if (!is.null(postscript.file)) 
        postscript(postscript.file, paper = "letter")
    resids <- 1:clusterStack$k
    if (sort.by == "resid" || !exists("p0")) 
        resids <- unlist(lapply(clusterStack, function(i) i$resid), 
            use.names = F)
    else if (sort.by == "p.clust") 
        resids <- unlist(lapply(clusterStack, function(i) i$p.clust), 
            use.names = F)
    else if (sort.by == "e.val") 
        resids <- unlist(lapply(clusterStack, function(i) mean(i$e.val, 
            na.rm = T)), use.names = F)
    else if (sort.by == "e.norm") 
        resids <- unlist(lapply(clusterStack, function(i) log10(min(i$e.val, 
            na.rm = T)) + i$nrows), use.names = F)
    resids[is.na(resids)] <- 999
    resids[is.infinite(resids)] <- 999
    if (sort.by != "none") 
        sorted <- sort(resids, index.return = T)
    else {
        sorted <- list()
        sorted$ix <- 1:length(ks)
    }
    if (length(ks) < clusterStack$k) 
        n.best <- clusterStack$k
    if (n.best > clusterStack$k) 
        n.best <- clusterStack$k
    for (kk in sorted$ix[1:n.best]) {
        if (!(kk %in% ks)) 
            next
        cat.new("cluster:", kk, sort.by, ":", format(resids[kk], 
            digits = 4), "\n")
        clusterStack[[kk]]$k <- kk
        plotCluster.motif(clusterStack[[kk]], ratios, seqs = seqs, 
            postscript = NULL, min.e.val = min.e.val, short = short.names, 
            both = both.motifs, shade = shade, networks = networks, 
            iter = iter)
        if (is.null(postscript.file)) {
            tmp.menu <- menu(c(" next plot", " quit"))
            if (tmp.menu == 0 || tmp.menu == 2) 
                break
        }
    }
    if (!is.null(postscript.file)) 
        graphics.off()
}
plot.network <-
function (node.names, edges, node.coords = NULL, vizmap = NULL, 
    edge.sep.frac = 150, ...) 
{
    ow <- options("warn")
    options(warn = -1)
    if (is.null(vizmap)) 
        vizmap <- get.vizmap(node.names, edges)
    node.count <- length(node.names)
    if (is.null(node.coords)) {
        node.coords <- matrix(0, nrow = node.count, ncol = 2)
        node.coords[, 1] <- (runif(node.count) - 0.5) * 100
        node.coords[, 2] <- (runif(node.count) - 0.5) * 100
    }
    W <- diff(range(node.coords[, 1]))
    L <- diff(range(node.coords[, 2]))
    WL <- max(c(W, L))
    edge.sep = WL/edge.sep.frac
    plot(node.coords, type = "n", ...)
    n.edges <- 0
    if (!is.null(edges)) 
        n.edges <- nrow(edges)
    for (i in 1:n.edges) {
        if (n.edges == 0) 
            break
        node1 <- edges[i, 1]
        node2 <- edges[i, 2]
        if (node1 == node2) 
            next
        same.edges <- which((edges[, 1] == node1 & edges[, 2] == 
            node2) | (edges[, 1] == node2 & edges[, 2] == node1))
        edge <- c(edges[i, 1], edges[i, 2])
        if (length(same.edges) > 1) {
            ll <- length(same.edges)
            for (e in 1:ll) {
                if (all(edges[same.edges[e], ] == edge)) {
                  ind <- same.edges[e]
                  off <- norm.unit.vector(node.coords[edge, 1], 
                    node.coords[edge, 2])
                  xs <- node.coords[edge, 1] + (edge.sep * off[1] * 
                    (as.integer(ll/2) - e))
                  ys <- node.coords[edge, 2] + (edge.sep * off[2] * 
                    (as.integer(ll/2) - e))
                  if (is.null(vizmap$arrow) || vizmap$arrow[ind] <= 
                    0) 
                    lines(xs, ys, col = vizmap$ecol[ind], lty = vizmap$etyp[ind], 
                      lwd = vizmap$ewid[ind])
                  else arrows(xs[1], ys[1], xs[2], ys[2], code = vizmap$arrow[ind], 
                    length = vizmap$asize[ind], col = vizmap$ecol[ind], 
                    lty = vizmap$etyp[ind], lwd = vizmap$ewid[ind])
                }
            }
        }
        else {
            xs <- node.coords[edge, 1]
            ys <- node.coords[edge, 2]
            if (is.null(vizmap$arrow) || vizmap$arrow[i] <= 0) 
                lines(xs, ys, col = vizmap$ecol[i], lty = vizmap$etyp[i], 
                  lwd = vizmap$ewid[i])
            else arrows(xs[1], ys[1], xs[2], ys[2], code = vizmap$arrow[i], 
                length = vizmap$asize[i], col = vizmap$ecol[i], 
                lty = vizmap$etyp[i], lwd = vizmap$ewid[i])
        }
    }
    if (!is.null(vizmap$nshape) && is.matrix(vizmap$nshape)) {
        for (i in 1:nrow(node.coords)) {
            ii <- i
            if (length(vizmap$nshape) == 1) 
                ii <- 1
            polygon(vizmap$nshape[ii] + vizmap$nsize * matrix(node.coords[1, 
                ], ncol = 2, nrow = nrow(vizmap$nshape[ii]), 
                byrow = T), col = vizmap$nfill, border = vizmap$ncol, 
                lty = vizmap$ntyp)
        }
    }
    else {
        points(node.coords, pch = vizmap$nshape, cex = vizmap$nsize, 
            col = vizmap$ncol, bg = vizmap$nfill, ...)
    }
    text.coords <- node.coords
    text.coords[, 1] <- text.coords[, 1] + WL * vizmap$loff.x
    text.coords[, 2] <- text.coords[, 2] + WL * vizmap$loff.y
    text(text.coords, labels = node.names, xpd = NA, col = vizmap$lcol, 
        font = vizmap$lfont, cex = vizmap$lsize)
    options(ow)
    invisible(node.coords)
}
print.new <-
function (x, log = get.global("out.log"), ...) 
{
    print(x, ...)
    if (!is.null(log) && log != "") 
        try(capture.output(print(x, ...), file = log, append = T))
}
print.params <-
function (params) 
{
    cat.new("STARTING PARAMETERS:\n\n")
    for (name in sort(names(params))) {
        cat.new(sprintf("%25s", name), "=\t")
        if (class(params[[name]]) == "function") {
            next
        }
        else if (length(params[[name]]) <= 1) {
            cat.new(params[[name]], "\n")
        }
        else {
            if (is.null(names(params[[name]]))) 
                cat.new(which(!is.na(params[[name]])), ": ")
            else cat.new(names(params[[name]])[which(!is.na(params[[name]]))], 
                ": ")
            cat.new(params[[name]][!is.na(params[[name]])], "\n")
        }
    }
}
prune.clusterStack <-
function (clusterStack, cluster.row.floor) 
{
    cat.new("Pruning cluster stack: old size=", clusterStack$k)
    tmp.clusts <- list()
    for (i in 1:clusterStack$k) {
        if (is.null(clusterStack[[i]]) || is.null(clusterStack[[i]]$rows) || 
            clusterStack[[i]]$nrows <= cluster.row.floor || length(clusterStack[[i]]$rows) < 
            cluster.row.floor) 
            next
        tmp.clusts[[length(tmp.clusts) + 1]] <- clusterStack[[i]]
        cat.new(clusterStack[[i]]$k, "->", length(tmp.clusts), 
            "\n")
        tmp.clusts[[length(tmp.clusts)]]$k <- length(tmp.clusts)
    }
    tmp.clusts$k <- length(tmp.clusts)
    cat.new("\tnew size=", tmp.clusts$k, "\n")
    invisible(tmp.clusts)
}
readExpression <-
function () 
{
    cat("READING RATIOS:", ratios.fname, "\n")
    data <- read.table(gzfile(ratios.fname), header = T, as.is = T)
    ratios <- as.matrix(data[, 2:(length(data) - 1)])
    rownames(ratios) <- toupper(as.character(data$GENE))
    cat("READING LAMBDAS:", lambdas.fname, "\n")
    data <- read.table(gzfile(lambdas.fname), header = T, as.is = T)
    lambdas <- as.matrix(data[, 2:(length(data) - 1)])
    rownames(lambdas) <- toupper(as.character(data$GENE))
    colnames(ratios) <- gsub("\\.", "-", colnames(ratios), extended = F)
    colnames(lambdas) <- gsub("\\.", "-", colnames(lambdas), 
        extended = F)
    test.data <- read.table(gzfile(paste(data.dir, "/marray/final_uv_repair_genes_not_lam_filtered.txt.gz", 
        sep = "")), header = TRUE)
    start.col <- 3
    end.col <- 7
    ratios.olig <- as.matrix(test.data[, start.col:end.col])
    rownames(ratios.olig) <- toupper(as.character(test.data$DESCRIPTION))
    start.col <- 8
    end.col <- 12
    lambdas.olig <- as.matrix(test.data[, start.col:end.col])
    rownames(lambdas.olig) <- toupper(as.character(test.data$DESCRIPTION))
    gene.ids <- row.names(ratios.olig)
    gene.names <- as.character(test.data$GENE.NAME)
    names(gene.names) <- toupper(as.character(test.data$DESCRIPTION))
    rm(ratios.olig, lambdas.olig)
    if (FALSE) {
        cat("circ\n")
        test.data <- read.table(gzfile("data/marray/merge_unscaled_all_trimmed.clone.gz"), 
            header = TRUE, row.names = 2)
        start.col <- 2
        end.col <- 26
        ratios.circ <- as.matrix(test.data[, start.col:end.col])
        rownames(ratios.circ) <- toupper(rownames(ratios.circ))
        start.col <- 27
        end.col <- 51
        lambdas.circ <- as.matrix(test.data[, start.col:end.col])
        rownames(lambdas.circ) <- toupper(rownames(lambdas.circ))
        cat("iron\n")
        test.data <- read.table(gzfile("data/marray/halo_metal_response.mrna.gz"), 
            header = TRUE, row.names = 1)
        numCol <- (dim(test.data)[2] - 1)/2
        start.col <- 2
        end.col <- numCol + 1
        ratios.metal <- as.matrix(test.data[, start.col:end.col])
        rownames(ratios.metal) <- toupper(rownames(ratios.metal))
        start.col <- end.col + 1
        end.col <- start.col + numCol - 1
        lambdas.metal <- as.matrix(test.data[, start.col:end.col])
        rownames(lambdas.metal) <- toupper(rownames(lambdas.metal))
        cat("iron2\n")
        test.data <- read.table(gzfile("data/marray/halo_FeSO4_ts_06252004.mrna.gz"), 
            header = TRUE, row.names = 2)
        numCol <- (dim(test.data)[2] - 1)/2
        start.col <- 2
        end.col <- numCol + 1
        ratios.metal2 <- as.matrix(test.data[, start.col:end.col])
        rownames(ratios.metal2) <- toupper(rownames(ratios.metal2))
        start.col <- end.col + 1
        end.col <- start.col + numCol - 1
        lambdas.metal2 <- as.matrix(test.data[, start.col:end.col])
        rownames(lambdas.metal2) <- toupper(rownames(lambdas.metal2))
        cat("gene pert\n")
        test.data <- read.table(gzfile("data/marray/all_HO-LO-L-D_nobat_afsq2_1phr1.gz"), 
            header = TRUE, row.names = 1)
        numCol <- (dim(test.data)[2] - 1)/2
        start.col <- 2
        end.col <- numCol + 1
        ratios.gene <- as.matrix(test.data[, start.col:end.col])
        rownames(ratios.gene) <- toupper(rownames(ratios.gene))
        start.col <- end.col + 1
        end.col <- start.col + numCol - 1
        lambdas.gene <- as.matrix(test.data[, start.col:end.col])
        rownames(lambdas.gene) <- toupper(rownames(lambdas.gene))
        cat("gamma\n")
        test.data <- read.table(gzfile("data/marray/halo_gamma_all_20040504_215646.mrna.gz"), 
            header = TRUE, row.names = 1)
        numCol <- (dim(test.data)[2] - 1)/2
        start.col <- 2
        end.col <- numCol + 1
        ratios.gama <- as.matrix(test.data[, start.col:end.col])
        rownames(ratios.gama) <- toupper(rownames(ratios.gama))
        start.col <- end.col + 1
        end.col <- start.col + numCol - 1
        lambdas.gama <- as.matrix(test.data[, start.col:end.col])
        rownames(lambdas.gama) <- toupper(rownames(lambdas.gama))
        ratios <- cbind(ratios.olig, ratios.circ[gene.ids, ], 
            ratios.metal[gene.ids, ], ratios.gene[gene.ids, ], 
            ratios.gama[gene.ids, ], ratios.metal2[gene.ids, 
                ])
        lambdas <- cbind(lambdas.olig, lambdas.circ[gene.ids, 
            ], lambdas.metal[gene.ids, ], lambdas.gene[gene.ids, 
            ], lambdas.gama[gene.ids, ], lambdas.metal2[gene.ids, 
            ])
        tsStartStops <- c(2, 3, 4, 5, 6, 30, 116, 120, 121, 130, 
            131, 141)
        colMap <- makeColMap(tsStartStops, dim(ratios)[2], mode = "mixed")
        colMap[[132]]$del.t <- 1/60
        colMap[[133]]$del.t <- 5/60
        colMap[[134]]$del.t <- 5/60
        colMap[[135]]$del.t <- 5/60
        colMap[[136]]$del.t <- 5/60
        colMap[[137]]$del.t <- 5/60
        colMap[[138]]$del.t <- 15/60
        colMap[[139]]$del.t <- 40/60
        colMap[[140]]$del.t <- 80/60
        colMap[[141]]$del.t <- 160/60
        colMap[[3]]$del.t <- 30/60
        colMap[[5]]$del.t <- 30/60
        colMap[[117]]$del.t <- 20/60
        colMap[[118]]$del.t <- 30/60
        colMap[[119]]$del.t <- 180/60
        colMap[[120]]$del.t <- 720/60
        colMap[[122]]$del.t <- 10/60
        colMap[[123]]$del.t <- 10/60
        colMap[[124]]$del.t <- 10/60
        colMap[[125]]$del.t <- 10/60
        colMap[[126]]$del.t <- 20/60
        colMap[[127]]$del.t <- 60/60
        colMap[[128]]$del.t <- 120/60
        colMap[[129]]$del.t <- 240/60
        colMap[[130]]$del.t <- 480/60
    }
    cat("DIMS = ", dim(ratios), dim(lambdas), "\n")
    knockIn <- list()
    knockIn.names <- character()
    for (gene in gene.ids) {
        if (length(grep(paste("^", gene.names[gene], "__", sep = ""), 
            colnames(ratios), value = TRUE, perl = T)) > 0) {
            knockIn[[gene]] <- grep(paste("^", gene.names[gene], 
                "__", sep = ""), colnames(ratios), value = TRUE, 
                perl = T)
            cat("knock out for ", gene, gene.names[gene], knockIn[[gene]], 
                "\n")
            knockIn.names <- c(knockIn.names, gene)
        }
    }
    knockIn.index <- gene.ids[which(rownames(ratios) %in% knockIn.names)]
    n.genes.tot <- nrow(ratios)
    n.conds <- ncol(ratios)
    rm(start.col, end.col, test.data)
    rm(ratios.circ, ratios.olig, lambdas.olig, lambdas.circ, 
        lambdas.gene, lambdas.gama, ratios.gene, ratios.gama)
    rm(ratios.metal2, lambdas.metal2)
    rm(ratios.metal, lambdas.metal)
    cat("Total number of genes: ", n.genes.tot, ", and conditions: ", 
        n.conds, "\n")
    invisible(list(ratios = ratios, lambdas = lambdas, knockIn = knockIn, 
        knockIn.names = knockIn.names))
}
residOneClust <-
function (ratios, knockouts = FALSE, knockout.vec = NA, varNorm = FALSE) 
{
    n.rows <- nrow(ratios)
    n.cols <- ncol(ratios)
    d.rows <- apply(ratios, 1, mean, na.rm = TRUE)
    d.cols <- apply(ratios, 2, mean, na.rm = TRUE)
    d.all <- mean(d.rows, na.rm = TRUE)
    rij <- ratios + d.all
    for (i in 1:n.rows) rij[i, ] <- rij[i, ] - d.cols
    for (j in 1:n.cols) rij[, j] <- rij[, j] - d.rows
    if (!knockouts) {
        average.r <- mean(abs(rij), na.rm = TRUE)
        if (varNorm) {
            row.var <- mean(apply(ratios, 1, var, use = "pairwise.complete.obs"), 
                na.rm = T)
            if (row.var > maxRowVar) 
                row.var <- maxRowVar
            average.r <- average.r/row.var
        }
        return(average.r)
    }
    else {
        row.sums <- numeric(n.rows)
        for (i in 1:n.rows) {
            if (knockout.vec[[i]][1] > 0) {
                omit <- knockout.vec[[i]]
                row.sums[i] <- mean(abs(rij[i, -omit]), na.rm = TRUE)
            }
            else {
                row.sums[i] <- mean(abs(rij[i, ]), na.rm = TRUE)
            }
        }
        average.r <- mean(row.sums)
        if (varNorm) {
            row.var <- mean(apply(ratios, 1, var, use = "pairwise.complete.obs"), 
                na.rm = T)
            if (row.var > maxRowVar) 
                row.var <- maxRowVar
            average.r <- average.r/row.var
        }
        return(average.r)
    }
}
residOneClust.by.col <-
function (ratios, varNorm = TRUE) 
{
    n.rows <- nrow(ratios)
    n.cols <- ncol(ratios)
    d.rows <- apply(ratios, 1, mean, na.rm = TRUE)
    d.cols <- apply(ratios, 2, mean, na.rm = TRUE)
    d.all <- mean(d.rows, na.rm = TRUE)
    rij <- ratios + d.all
    for (i in 1:n.rows) rij[i, ] <- rij[i, ] - d.cols
    for (j in 1:n.cols) rij[, j] <- rij[, j] - d.rows
    average.r <- apply(abs(rij), 2, mean, na.rm = TRUE)
    if (varNorm) {
        row.var <- mean(apply(ratios, 1, var, use = "pairwise.complete.obs"), 
            na.rm = T)
        if (row.var > maxRowVar) 
            row.var <- maxRowVar
        average.r <- average.r/row.var
    }
    return(average.r)
}
residOneClust.noCap <-
function (ratios, knockouts = FALSE, knockout.vec = NA, varNorm = FALSE) 
{
    n.rows <- nrow(ratios)
    n.cols <- ncol(ratios)
    d.rows <- apply(ratios, 1, mean, na.rm = TRUE)
    d.cols <- apply(ratios, 2, mean, na.rm = TRUE)
    d.all <- mean(d.rows, na.rm = TRUE)
    rij <- ratios + d.all
    for (i in 1:n.rows) rij[i, ] <- rij[i, ] - d.cols
    for (j in 1:n.cols) rij[, j] <- rij[, j] - d.rows
    if (!knockouts) {
        average.r <- mean(abs(rij), na.rm = TRUE)
        if (varNorm) {
            row.var <- mean(apply(ratios, 1, var, use = "pairwise.complete.obs"))
            average.r <- average.r/row.var
        }
        return(average.r)
    }
    else {
        row.sums <- numeric(n.rows)
        for (i in 1:n.rows) {
            if (knockout.vec[[i]][1] > 0) {
                omit <- knockout.vec[[i]]
                row.sums[i] <- mean(abs(rij[i, -omit]), na.rm = TRUE)
            }
            else {
                row.sums[i] <- mean(abs(rij[i, ]), na.rm = TRUE)
            }
        }
        average.r <- mean(row.sums)
        if (varNorm) {
            row.var <- mean(apply(ratios, 1, var, use = "pairwise.complete.obs"))
            average.r <- average.r/row.var
        }
        return(average.r)
    }
}
rev.comp <-
function (seqs) 
{
    out <- gsub("A", "t", seqs)
    out <- gsub("T", "a", out)
    out <- gsub("C", "g", out)
    out <- toupper(gsub("G", "c", out))
    names(out) <- names(seqs)
    out <- sapply(out, function(i) paste(rev(strsplit(i, "")[[1]]), 
        sep = "", collapse = ""))
    names(out) <- names(seqs)
    return(out)
}
rowVar <-
function (ratios, clust, mode = "all") 
{
    row.var <- numeric()
    if (mode == "all") {
        for (k in 1:clust$k) {
            row.var[k] <- mean(apply(ratios[clust[[k]]$rows, 
                clust[[k]]$cols], 1, var, na.rm = T))
        }
    }
    else if (mode > 0 && mode <= clust$k) {
        k <- mode
        row.var <- mean(apply(ratios[clust[[k]]$rows, clust[[k]]$cols], 
            1, var, na.rm = T))
    }
    else {
        stop("impropper mode specified... should be all or 1:numClusts\n")
    }
    return(row.var)
}
run <-
function () 
{
    source("cMonkey.R")
}
runMast <-
function (memeOut, genes, seqs, fname = "meme.tmp.fst", memeOutFname = "meme.tmp.out", 
    bgseqs = NULL, bg.list = NULL, bgfname = NULL, e.value.cutoff = 99999, 
    p.value.cutoff = 0.5, motif.e.value.cutoff = 99999, unlink = T, 
    verbose = F) 
{
    fname <- tempfile(fname)
    if (!is.null(bgfname)) 
        bgfname <- tempfile(bgfname)
    memeOutFname <- tempfile(memeOutFname)
    cat(memeOut, sep = "\n", file = memeOutFname)
    lengths <- mkTempMemeFiles(genes, seqs, fname = fname, bgseqs = bgseqs, 
        bg.list = bg.list, bgfname = bgfname)
    cmd <- paste(mast.cmd, memeOutFname, "-d", fname, "-bfile", 
        bgfname, "-nostatus -stdout -text -seqp -remcorr -ev ", 
        e.value.cutoff, " -mev ", motif.e.value.cutoff, " -mt ", 
        p.value.cutoff, " -brief")
    if (verbose) 
        cat.new(cmd, "\n")
    output <- system(cmd, intern = TRUE)
    if (unlink) 
        my.unlink(c(fname, memeOutFname))
    output
}
runMeme <-
function (sgenes, seqs, fname = "meme.tmp.fst", bgseqs = NULL, 
    bgfname = NULL, bg.list = NULL, nmotif = 1, e.value.cutoff = 99999, 
    unlink = T, minw = 6, maxw = 17, model = "zoops", verbose = F) 
{
    fname <- tempfile(fname)
    if (!is.null(bgfname)) 
        bgfname <- tempfile(bgfname)
    lengths <- mkTempMemeFiles(sgenes, seqs, fname = fname, bgseqs = bgseqs, 
        bg.list = bg.list, bgfname = bgfname)
    cmd <- paste(meme.cmd, fname, "-bfile", bgfname, "-nostatus -dna -revcomp -text -nmotifs", 
        nmotif, "-mod", model, "-minw", minw, "-maxw", maxw, 
        "-time 1000 -evt", e.value.cutoff, "-maxsize", lengths)
    if (verbose) 
        cat.new(cmd, "\n")
    output <- system(cmd, intern = TRUE, ignore.stderr = TRUE)
    if (unlink) 
        my.unlink(c(fname))
    return(output)
}
space.pad <-
function (lines, length) 
{
    for (i in 1:length(lines)) {
        nc <- nchar(lines[i])
        if (nc >= length) 
            next
        lines[i] <- paste(c(lines[i], rep(" ", length - nc)), 
            sep = "", collapse = "")
    }
    return(lines)
}
start.up <-
function () 
{
    initialize.biclust(T)
}
subnet.p.value <-
function (genes, network, gene.ids, edge.centric = T, n.genes = NULL, 
    min.edges = 3, adj.list = NULL) 
{
    p.value <- 1
    if (edge.centric) {
        nr <- length(gene.ids)
        total.edges <- nrow(network)/2
        total.possible.edges <- (nr * nr - nr)/2
        ng <- length(genes)
        tmp <- unlist(adj.list[genes])
        if (is.na(total.edges) || total.edges < min.edges) 
            return(1.01)
        sub.edges <- 0
        if (!is.null(tmp) && length(tmp) > 0) 
            sub.edges <- sum(tmp %in% genes)
        if (sub.edges < min.edges) 
            return(2)
        sub.possible.edges <- (ng * ng - ng)/2
        p.value <- phyper(sub.edges, total.edges, total.possible.edges - 
            total.edges, sub.possible.edges, lower.tail = F)
    }
    else {
        stop("NOT SUPPORTED\n")
    }
    if (p.value <= 1e-30) 
        p.value <- 1e-30
    return(p.value)
}
subnet.p.values <-
function (genes, networks, gene.ids, edge.centric = T, net.q0 = NA, 
    force = F) 
{
    p.values <- numeric()
    for (i in names(networks)) {
        if (i == "n.edges" || i == "n.genes" || i == "n.conds" || 
            i == "cond.sims" || i == "sums" || i == "adj.lists") 
            next
        if (!force && !is.na(net.q0) && (is.na(net.q0[i]) || 
            net.q0[i] <= 0)) {
            p.values[i] <- 2
            next
        }
        p.values[i] <- subnet.p.value(genes, networks[[i]], gene.ids, 
            edge.centric, n.genes = networks$n.genes, adj.list = networks$adj.lists[[i]])
    }
    p.values
}
subSeq <-
function (seqs, len = 100, direction = "back") 
{
    if (nchar(seqs[1]) < len) 
        return(seqs)
    new.seqs <- character()
    if (direction == "back" || direction == "b") {
        len <- len - 1
        for (i in 1:length(seqs)) new.seqs[i] <- substr(seqs[i], 
            nchar(seqs[i]) - len, nchar(seqs[i]))
    }
    else if (direction == "forward" || direction == "f") {
        for (i in 1:length(seqs)) new.seqs[i] <- substr(seqs[i], 
            1, len)
    }
    names(new.seqs) <- names(seqs)
    return(new.seqs)
}
symmetrize.edge.list <-
function (edges, remove.diag = T) 
{
    out <- uniquify.edge.list(edges)
    out <- rbind(out, cbind(out[, 2], out[, 1]))
    invisible(out)
}
symmetrize.net.matrix <-
function (mat) 
{
    if (sum(mat) == 0) 
        return(mat)
    tmp <- which(mat != 0, arr = T)
    tmp2 <- cbind(tmp[, "col"], tmp[, "row"])
    out <- mat
    out[tmp2] <- mat[tmp]
    return(out)
}
test.opt.all.clusters <-
function (clusterStack, ks = 1:get.global("kmax"), debug = F, 
    mean.count = NA, ...) 
{
    if (length(clusterStack) == 0) 
        clusterStack$k <- 0
    else for (k in 1:clusterStack$k) clusterStack[[k]]$k <- k
    if (!exists("upstream.bg") || is.null(upstream.bg) || upstream.bg$order != 
        bg.order) {
        upstream.bg <<- mkBgFile(upstream, order = bg.order, 
            bgfname = paste(data.dir, "/upstream.bg", sep = ""))
    }
    cl.size <- 1
    samp <- ks
    if (clusterStack$k > 0) 
        samp <- samp[!samp %in% 1:clusterStack$k]
    if (!exists("out.logs")) 
        out.logs <- list()
    sampled.genes <- character()
    if (length(clusterStack) > 1) 
        sampled.genes <- unique(sapply(clusterStack[1:clusterStack$k], 
            "[[", "rows"))
    for (ind in seq(min(samp), max(samp), cl.size)) {
        inds <- ind:(ind + cl.size - 1)
        clusterStack <- addClusters.semirnd(inds, clusterStack, 
            ratios, no.overlap = T, clust.size = semirnd.clust.size, 
            dont.allow = sampled.genes)
        sampled.genes <- unique(c(sampled.genes, unlist(lapply(clusterStack[1:clusterStack$k], 
            "[[", "rows"))))
        counts <- NULL
        if (!is.na(mean.count) && clusterStack$k > 1) {
            if (min(inds) > 1) {
                counts <- genes.in.how.many.clusters(clusterStack, 
                  gene.ids, up.to = (min(inds) - 1))
                cat.new(sum(counts == 0), "genes not in any clusters... yet!\n")
                if (max(counts) > 2) {
                  tmp.hist <- hist(counts, plot = F, breaks = max(counts))
                  tmp.x <- tmp.hist$counts + 1
                  names(tmp.x) <- as.character(tmp.hist$breaks[1:length(tmp.x)])
                  cat.new("Histogram of number of genes in each cluster:\n")
                  print.new(tmp.x)
                  rm(tmp.hist, tmp.x)
                  gc()
                }
            }
        }
        cat.new(date(), "\n")
        cat.new("Optimizing these clusters:", inds, "\n")
        ll <- lapply(clusterStack[inds], test.opt.one.cluster, 
            plot = F, plot.cluster = F, log = T, counts = counts, 
            mean.count = mean.count, ...)
        for (ii in 1:length(ll)) {
            i <- inds[ii]
            new.clust <- ll[[ii]]$cluster
            log <- ll[[ii]]$log
            if (class(new.clust) != "list") {
                cat.new(ind, ind + ii - 1, i, "Error: ", new.clust, 
                  "\n")
                next
            }
            if (is.null(new.clust) || is.na(i)) 
                next
            out.logs[[i]] <- log
            cat.new("OLD/NEW SCORES FOR CLUSTER:", ind, ind + 
                ii - 1, i, ":   ", clusterStack[[i]]$nrows, new.clust$nrows, 
                "\t", clusterStack[[i]]$resid, new.clust$resid, 
                "\t", new.clust$p.clust, "\t", min(new.clust$motif.out$e.value, 
                  na.rm = T), "\n")
            new.clust$motif.out$diagrams <- new.clust$motif.out$diagrams[new.clust$rows]
            new.clust$motif.out$p.values <- new.clust$motif.out$p.values[new.clust$rows]
            if (new.clust$nrows < cluster.row.floor) {
                new.clust$nrows <- 0
                new.clust$rows <- character()
            }
            clusterStack[[i]] <- new.clust
            f1 <- paste(output.dir, "/clust_", sprintf("%03d", 
                as.integer(i)), ".ps", sep = "")
            f1.pdf <- sub(".ps$", ".pdf", f1, perl = T)
            if (file.exists(f1)) {
                system(paste(ps2pdf.bin, f1, f1.pdf))
                unlink(f1)
            }
            f2 <- paste(output.dir, "/cstats_", sprintf("%03d", 
                as.integer(i)), ".ps", sep = "")
            f2.pdf <- sub(".ps$", ".pdf", f2, perl = T)
            if (file.exists(f2)) {
                system(paste(ps2pdf.bin, f2, f2.pdf))
                unlink(f2)
            }
        }
        curr.state.file <- paste(state.file, ".", sprintf("%03d", 
            as.integer(i)), ".RData", sep = "")
        save(clusterStack, out.logs, date.biclust.run, biclust.version, 
            file = curr.state.file, compress = T)
    }
    clusterStack$k <- max(ks)
    invisible(list(clusts = clusterStack, logs = out.logs))
}
test.opt.one.cluster <-
function (cluster, bad.move.temp = seq(0.18, 0.1, length = 20), 
    max.moves = 5, iters = 1:20, rows.only = F, cols.only = F, 
    add.max = 3, remove.max = 3, max.rows = 15, weight.fixed = 1, 
    counts = NULL, mean.count = 2, sample = T, update.params = T, 
    log.it = T, plot = T, plot.postscript = NULL, plot.cluster = F, 
    plot.cluster.postscript = NULL) 
{
    log <- list()
    k <- cluster$k
    cat.new(0, k, "RESID:", cluster$resid, "\tPCLUST:", cluster$p.clust, 
        "\n")
    CLUSTER <<- cluster
    use.nets <- NULL
    use.nets <- NULL
    if (!is.null(cluster$net.p.old) && any(net.q0 > 0)) 
        use.nets <- names(which(net.q0 > 0))[which(names(which(net.q0 > 
            0)) %in% names(cluster$net.p.old))]
    if (!is.null(use.nets)) {
        cat.new(0, k, "SUBNET PVALS:\n")
        print.new(cluster$net.p.old[use.nets])
    }
    iters <- sort(iters)
    if (length(bad.move.temp) < length(iters)) 
        bad.move.temp <- rep(bad.move.temp, length = length(iters))
    if (log.it) 
        log[["0"]] <- list(nrows = cluster$nrows, ncols = cluster$ncols, 
            resid = cluster$resid, p.clust = cluster$p.clust, 
            net.pvals = cluster$net.p.old, log.lik = NA)
    dir.name <- output.dir
    if (plot || plot.cluster && !file.exists(dir.name)) 
        try(dir.create(dir.name, showWarnings = F))
    dir.name <- paste(ps2pdf.cmd, output.dir)
    if (plot) {
        plot.postscript <- paste(dir.name, "/cstats_", sprintf("%03d", 
            as.integer(k)), ".pdf", sep = "")
    }
    if (plot.cluster) {
        plot.cluster.postscript <- paste(dir.name, "/clust_", 
            sprintf("%03d", as.integer(k)), ".pdf", sep = "")
    }
    plot.dev <- plot.cluster.dev <- NULL
    changed <- TRUE
    old.n.motifs <- -1
    for (iter in iters) {
        gc()
        old.n.motifs <- n.motifs
        old.changed <- changed
        if (update.params) {
            try(detach("local.iter.params"), silent = T)
            local.iter.params <- get.all.iter.based.params(iter, 
                params, verbose = F)
            attach(local.iter.params, pos = -1)
        }
        cat.new(iter, k, "r0:", r0, " p0:", p0, " q0:", q0, " n.motifs:", 
            n.motifs, " TEMP:", bad.move.temp[iter], "\n")
        if (plot.cluster) {
            if (is.null(plot.cluster.dev)) {
                if (is.null(plot.cluster.postscript)) {
                  plot.cluster.dev <- 3
                }
                else {
                  postscript(plot.cluster.postscript, paper = "letter")
                  plot.cluster.dev <- dev.cur()
                }
            }
            if (old.changed || old.n.motifs != n.motifs) {
                dev.set(plot.cluster.dev)
                try(plotCluster.motif(cluster, iter = iter), 
                  silent = T)
            }
        }
        if (is.null(plot.dev)) {
            if (is.null(plot.postscript)) {
                plot.dev <- 2
            }
            else {
                postscript(plot.postscript, paper = "letter")
                plot.dev <- dev.cur()
            }
        }
        dev.set(plot.dev)
        par(mfrow = c(3, 3))
        changed <- FALSE
        if (!cols.only) {
            if (old.changed || old.n.motifs != n.motifs) {
                rows <- cluster$rows
                not.rows <- gene.ids[!(gene.ids %in% rows)]
                find.best.weights <- T
                if (length(rows) > max.rows) 
                  find.best.weights <- F
                row.gains <- get.row.gains.for.cluster(cluster, 
                  plot = plot, bad.move.temp = bad.move.temp[iter], 
                  response = T, find.best.weights = find.best.weights, 
                  weight.fixed = weight.fixed)
            }
            resp <- row.gains
            log.lik <- get.cluster.likelihood(cluster, resp)
            if (sample) {
                sample.prob <- exp((resp - 1)/bad.move.temp[iter])
                if (!is.null(counts)) {
                  count.prob <- ppois(counts, mean.count, lower = F)/ppois(mean.count, 
                    mean.count, lower = F)
                  count.prob[rows] <- ppois(counts[rows], mean.count, 
                    lower = T)/ppois(mean.count, mean.count, 
                    lower = T)
                  sample.prob <- sample.prob * count.prob
                }
                if (plot && (old.changed || old.n.motifs != n.motifs)) {
                  plot(row.gains, sample.prob, main = paste("iter=", 
                    iter, "; temp=", bad.move.temp[iter]))
                  points(row.gains[rows], sample.prob[rows], 
                    pch = 20, col = "red")
                }
                samples <- runif(length(sample.prob))
                allowed.moves <- sample.prob > samples
            }
            else {
                allowed.moves <- (row.gains > 0)
            }
            if (sum(allowed.moves, na.rm = T) > 0) {
                hist(row.gains[allowed.moves], breaks = 50)
                if (sum(allowed.moves, na.rm = T) > max.moves) {
                  tmp <- sort(resp[allowed.moves], dec = T)[1:max.moves]
                  allowed.moves <- gene.ids %in% names(tmp)
                  names(allowed.moves) <- gene.ids
                }
                to.move <- names(which(allowed.moves == TRUE))
                to.remove <- to.move[to.move %in% rows]
                if (length(to.remove) > remove.max) 
                  to.remove <- to.remove[1:remove.max]
                to.add <- to.move[!(to.move %in% rows)]
                if (length(to.add) > add.max) 
                  to.add <- to.add[1:add.max]
                if (length(to.add) > 0) {
                  cat.new(iter, k, "ROWS TO ADD:", to.add, row.gains[to.add], 
                    "\n")
                  cluster$rows <- unique(c(cluster$rows, to.add))
                  if (!is.null(counts)) 
                    counts[to.add] <- counts[to.add] + 1
                }
                if (length(to.remove) > 0) {
                  cat.new(iter, k, "ROWS TO REMOVE:", to.remove, 
                    row.gains[to.remove], "\n")
                  cluster$rows <- cluster$rows[!(cluster$rows %in% 
                    to.remove)]
                  if (!is.null(counts)) 
                    counts[to.remove] <- counts[to.remove] - 
                      1
                }
                cluster$nrows <- length(cluster$rows)
                if (cluster$nrows < cluster.row.floor) {
                  cat.new("Cluster has shrunk too small... breaking out.\n")
                  if (log.it) {
                    for (tmp.iter in (iter + 1):max(iters)) log[[as.character(tmp.iter)]] <- list(nrows = cluster$nrows, 
                      ncols = cluster$ncols, resid = cluster$resid, 
                      p.clust = cluster$p.clust, net.pvals = cluster$net.p.old, 
                      log.lik = log.lik)
                  }
                  break
                }
                if (length(to.add) + length(to.remove) > 0) {
                  cluster$motif.out <- motif.one.cluster(cluster, 
                    n.motifs = n.motifs)
                  cluster$motif.out$mast.out <- NULL
                  cluster <- fillPClust(cluster, k = cluster$k, 
                    e.val.cutoff = e.val.cutoff, p.val.cutoff = p.val.cutoff)
                  cluster$net.p.old <- get.cluster.subnet.p.values(cluster, 
                    networks, gene.ids, edge.centric = T, net.q0 = net.q0, 
                    rows.only = T, force = F)
                  changed <- TRUE
                }
            }
        }
        if (!rows.only) {
            if (changed || old.changed || old.n.motifs != n.motifs) {
                cols <- cluster$cols
                not.cols <- all.conditions[!(all.conditions %in% 
                  cols)]
                col.gains <- get.col.gains.for.cluster(cluster, 
                  bad.move.temp = bad.move.temp[iter], plot = plot, 
                  response = T)
            }
            resp <- col.gains
            if (sample) {
                sample.prob <- exp((resp - 1)/bad.move.temp[iter])
                samples <- runif(length(sample.prob))
                allowed.moves <- sample.prob > samples
                if (plot && (old.changed || old.n.motifs != n.motifs)) {
                  plot(col.gains, sample.prob, main = paste("iter=", 
                    iter))
                  points(col.gains[cols], sample.prob[cols], 
                    pch = 20, col = "red")
                }
            }
            else {
                allowed.moves <- col.gains > 0
            }
            if (sum(allowed.moves, na.rm = T) > 0) {
                if (sum(allowed.moves, na.rm = T) > max.moves) {
                  tmp <- sort(resp[allowed.moves], dec = T)[1:max.moves]
                  allowed.moves <- all.conditions %in% names(tmp)
                  names(allowed.moves) <- all.conditions
                }
                to.move <- names(which(allowed.moves == TRUE))
                to.remove <- to.move[to.move %in% cols]
                to.add <- to.move[!(to.move %in% cols)]
                if (length(to.add) > 0) {
                  cat.new(iter, k, "COLS TO ADD:", to.add, col.gains[to.add], 
                    "\n")
                  cluster$cols <- unique(c(cluster$cols, to.add))
                }
                if (length(to.remove) > 0) {
                  cat.new(iter, k, "COLS TO REMOVE:", to.remove, 
                    col.gains[to.remove], "\n")
                  cluster$cols <- cluster$cols[!(cluster$cols %in% 
                    to.remove)]
                }
                cluster$ncols <- length(cluster$cols)
                if (length(to.add) + length(to.remove) > 0) 
                  changed <- TRUE
            }
        }
        if (changed || old.n.motifs != n.motifs) {
            cluster$resid <- residOneClust(ratios[cluster$rows, 
                cluster$cols], varNorm = TRUE)
            if (net.q0["cond.sims"] > 0) 
                cluster$net.p.old["cond.sims"] <- condition.sim.pvalue(cluster$cols, 
                  networks, all.conditions)
        }
        cat.new(iter, k, "NROWS:", cluster$nrows, "\tNCOLS:", 
            cluster$ncols)
        if (is.na(cluster$motif.out$e.value)) 
            cat.new("\tRESID:", cluster$resid, "\tPCLUST:", cluster$p.clust, 
                "\n")
        else cat.new("\tRESID:", cluster$resid, "\tPCLUST:", 
            cluster$p.clust, "\tEVAL:", min(cluster$motif.out$e.value, 
                na.rm = T), "\tLOGLIK:", log.lik, "\n")
        if (!is.null(use.nets)) {
            cat.new(iter, k, "SUBNET PVALS:\n")
            print.new(cluster$net.p.old[use.nets])
        }
        if (log.it) 
            log[[as.character(iter)]] <- list(nrows = cluster$nrows, 
                ncols = cluster$ncols, resid = cluster$resid, 
                p.clust = cluster$p.clust, e.values = cluster$motif.out$e.value, 
                net.pvals = cluster$net.p.old, log.lik = log.lik)
        if (cluster$nrows < cluster.row.floor) {
            cat.new("Cluster is invalid... breaking out.\n")
            if (log.it) {
                for (tmp.iter in (iter + 1):max(iters)) log[[as.character(tmp.iter)]] <- list(nrows = cluster$nrows, 
                  ncols = cluster$ncols, resid = cluster$resid, 
                  p.clust = cluster$p.clust, net.pvals = cluster$net.p.old, 
                  log.lik = log.lik)
            }
            break
        }
    }
    if (!is.null(plot.postscript)) {
        dev.set(plot.dev)
        par(mfrow = c(3, 3))
        try(plot(sapply(log, "[[", "resid"), type = "l", main = paste("Cluster", 
            k, "resids"), xlab = "iteration"))
        try(plot(sapply(log, "[[", "p.clust"), col = "red", type = "l", 
            main = paste("Cluster", k, "p.clust"), xlab = "iteration"), 
            silent = T)
        try(plot(sapply(log, "[[", "nrows"), col = "green", type = "l", 
            main = paste("Cluster", k, "nrows"), xlab = "iteration"), 
            silent = T)
        try(plot(sapply(log, "[[", "ncols"), col = "blue", type = "l", 
            main = paste("Cluster", k, "ncols"), xlab = "iteration"), 
            silent = T)
        try(plot(apply(sapply(log, "[[", "e.values"), min, na.rm = T), 
            col = "violet", type = "l", main = paste("Cluster", 
                k, "e.value"), xlab = "iteration"), silent = T)
        try(plot(sapply(log, "[[", "log.lik"), col = "red", type = "l", 
            main = paste("Cluster", k, "log likelihood"), xlab = "iteration"), 
            silent = T)
        if (is.null(log[[2]]$net.pvals)) 
            log[[2]]$net.pvals <- log[[3]]$net.pvals
        if (is.null(log[[1]]$net.pvals)) 
            log[[1]]$net.pvals <- log[[2]]$net.pvals
        nets <- log10(sapply(log, "[[", "net.pvals"))
        for (net.name in use.nets) try(plot(nets[net.name, ], 
            col = "orange", type = "l", main = paste("Cluster", 
                k, net.name, "log10-p-values"), xlab = "iteration"), 
            silent = T)
    }
    if (!is.null(plot.postscript)) 
        dev.off(plot.dev)
    if (!is.null(plot.cluster.postscript)) 
        dev.off(plot.cluster.dev)
    try(detach("local.iter.params"), silent = T)
    if (!log.it) 
        invisible(cluster)
    else invisible(list(cluster = cluster, log = log))
}
uniquify.edge.list <-
function (edges, remove.diag = T) 
{
    if (is.vector(edges)) {
        edges <- matrix(edges, nrow = 1)
        return(edges)
    }
    edges <- unique(edges)
    if (is.vector(edges)) {
        edges <- matrix(edges, nrow = 1)
        return(edges)
    }
    out <- matrix(edges[1, ], ncol = 2)
    self.edges <- numeric()
    if (remove.diag) 
        self.edges <- which(edges[, 1] == edges[, 2])
    for (i in 2:nrow(edges)) {
        if (i %in% self.edges) 
            next
        tmp1 <- out[, 1] %in% edges[i, ]
        if (sum(tmp1) == 0) {
            out <- rbind(out, edges[i, ])
            next
        }
        tmp2 <- out[tmp1, 2] %in% edges[i, ]
        if (sum(tmp2) == 0) 
            out <- rbind(out, edges[i, ])
    }
    invisible(out)
}
viewPssm <-
function (pssm, e.val = "rnd", mode = c("scaled", "raw"), mot.ind = 1, 
    separate.pages = F) 
{
    if (is.null(pssm)) 
        return()
    colMap <- c(3, 4, 7, 2)
    colLet <- c("A", "C", "G", "T")
    win.size <- dim(pssm)[1]
    old.pars <- par()
    if (!separate.pages) 
        par(mar = rep(0.5, 4), mgp = c(3, 1, 0) * 0.5)
    if (mode == "raw") {
        x.range <- c(0.5, win.size + 0.5)
        y.range <- c(0, 1)
        plot(x.range, y.range, type = "n")
        for (j in 1:win.size) {
            inds <- sort(pssm.sc[j, ], index = T)$ix
            for (i in 1:4) {
                if (i == 1) {
                  rect((j - 0.5), 0, (j + 0.5), pssm[j, ind], 
                    col = colMap[ind])
                  if (pssm[j, ind] > 0.05) 
                    text(j, 0 + pssm[j, ind]/2, colLet[ind])
                  prev.h <- pssm[j, ind]
                }
                else {
                  rect((j - 0.5), prev.h, (j + 0.5), (pssm[j, 
                    ind] + prev.h), col = colMap[ind])
                  if (pssm[j, ind] > 0.05) {
                    if (i == 2) 
                      text(j, prev.h + 0.5 * pssm[j, ind], colLet[ind], 
                        col = 8)
                    else text(j, prev.h + 0.5 * pssm[j, ind], 
                      colLet[ind])
                  }
                  prev.h <- prev.h + pssm[j, i]
                }
            }
        }
        tmp.tit <- paste("Raw PSSM #", mot.ind, ": E=", e.val, 
            sep = "")
        title(tmp.tit, cex.main = 0.7)
    }
    else {
        entr <- getEntropy(pssm)
        scale.e <- (2 - entr)/2
        scale.e[scale.e < 0.05] <- 0.05
        x.range <- c(0.5, win.size + 0.5)
        if (max(scale.e) > 0.7) 
            y.range <- c(0, 1)
        else y.range <- c(0, 0.7)
        plot(x.range, y.range, type = "n", tck = 0.01, cex.lab = 0.2, 
            cex.sub = 0.2, cex.axis = 0.2, axes = F)
        tmp.tit <- paste("Scaled PSSM #", mot.ind, ": E=", e.val, 
            sep = "")
        title(tmp.tit, cex.main = 0.7)
        pssm.sc <- scale.e * pssm
        for (j in 1:win.size) {
            inds <- sort(pssm.sc[j, ], index = T)$ix
            for (i in 1:4) {
                ind <- inds[i]
                if (i == 1) {
                  rect((j - 0.5), 0, (j + 0.5), pssm.sc[j, ind], 
                    col = colMap[ind])
                  if (pssm[j, ind] > 0.05) 
                    text(j, 0 + pssm.sc[j, ind]/2, colLet[ind])
                  prev.h <- pssm.sc[j, ind]
                }
                else {
                  rect((j - 0.5), prev.h, (j + 0.5), (pssm.sc[j, 
                    ind] + prev.h), col = colMap[ind])
                  if (pssm.sc[j, ind] > 0.05) {
                    if (i == 2) 
                      text(j, prev.h + 0.5 * pssm.sc[j, ind], 
                        colLet[ind], col = 8)
                    else text(j, prev.h + 0.5 * pssm.sc[j, ind], 
                      colLet[ind])
                  }
                  prev.h <- prev.h + pssm.sc[j, ind]
                }
            }
        }
        if (win.size < 10) 
            text(1:win.size, rep(-0.01, win.size), as.character(1:win.size), 
                cex = 0.7, adj = c(0.5, 1), xpd = NA)
        else text(seq(1, win.size, 2), rep(-0.01, win.size), 
            as.character(seq(1, win.size, 2)), cex = 0.7, adj = c(0.5, 
                1), xpd = NA)
    }
    lines(c(0, win.size + 0.5), c(0.05, 0.05), col = 1, lty = 2, 
        cex = 3, lwd = 3)
}
