#' S4 class containing the data to draw a regional association plot with the 
#'   'ld' style.
#' 
#' @slot ldPValues Data frame containing the p-values for the 'ld' style.
#' @slot indexVariant Data frame with one row, containing the index variant.
#' 
#' @seealso \code{\link{regionplot}}
styleLd <- setClass("styleLd", 
                    slots = c(ldPvalues = 'data.frame',
                              indexVariant = 'data.frame'),
                    prototype = c(ldPValues = NULL,
                                  indexVariant = NULL))

#' S4 class containing the data to draw a regional association plot with the 
#'   'signals' style.
#'   
#' @slot signalsPValues Data frame containing the p-values for the 'signals' 
#'   style.
#' @slot index_cleo UNDEFINED
#' 
#' @seealso \code{\link{regionplot}}
styleSignals <- setClass("styleSignals",
                         slots = c(signalsPValues = 'data.frame',
                                   index_cleo = 'data.frame'),
                         prototype = c(signalsPValues = NULL,
                                       index_cleo = NULL))

#' S4 class containing the data to draw a regional association plot with the 
#'   'signal' style.
#'   
#' @slot signalPValues Data frame with p-values for the 'signal' style.
#' 
#' @seealso \code{\link{regionplot}}
styleSignal <- setClass("styleSignal",
                        slots = c(signalPValues = 'data.frame'),
                        prototype = c(signalPValues = NULL))

#' S4 class containing the data to draw a regional association plot.
#'  
#' @slot queryForPValues Query used to get the p-values for a given GWAS 
#'   from a database.
#' @slot pValues Data frame containing the original p-values.
#' @slot xEntity Named list containing \code{entity}, \code{entity_label} and 
#'   \code{entity_where}.
#' @slot xRegion Named list containing \code{chrom}, \code{pos_start}, 
#'   \code{pos_end} and \code{label}. 
#'   \itemize{
#'    \item{\code{chrom} is the chromosome to be plotted.} 
#'    \item{\code{pos_start} and \code{pos_end} specify the start and end of 
#'     the region to draw the plot for.} 
#'    \item{\code{label} gives a label for printing that shows the chromosome 
#'     and region.}
#'    }
#' @slot indexSNPs Data frame with two columns for the position and p-values of 
#'   SNPs that should be highlighted in the output plot.
#' @slot plotTitle Title for a regional association plot.
#' @slot plotSubTitle Subtitle for a regional association plot.
#' @slot styleSignals Named list with the data needed to draw a regional 
#'   association plot with the 'signals' style. List contains:
#'   \itemize{
#'    \item{\code{signalsPvalues} is the p-values data frame for style 
#'      'signals'.}
#'    \item{\code{indexCleo} is the summary statistics for index variants, which
#'      is generated during CLEO finemapping.}
#'    } 
#' @slot styleSignal Data frame with the p-values for style 'signal'.
#' 
#' @seealso \code{\link{regionplot}}
regionplot <- setClass("regionplot",
                       slots = c(queryForPValues = "character",
                                 pValues = 'data.frame',
                                 xEntity = "list",
                                 xRegion = "list",
                                 indexSNPs = 'data.frame',
                                 plotTitle = "character",
                                 plotSubTitle = "character",
                                 styleSignals = 'styleSignals',
                                 styleSignal = 'data.frame',
                                 styleLd = 'styleLd'),
                       prototype = list(queryForPValues = NA_character_,
                                        pValues = NULL,
                                        xEntity = NULL,
                                        xRegion = NULL,
                                        indexSNPs = NULL,
                                        plotTitle = NA_character_,
                                        plotSubTitle = NA_character_,
                                        styleSignals = NULL,
                                        styleSignal = NULL,
                                        styleLd = NULL)
)

#' Regional association plot.
#' 
#' Draw a regional association plot.
#' 
#' @param analysis Character giving the key value for the GWAS analysis to draw 
#'   a regional association plot from.
#' @param entity Optional entity within eQTL to draw the regional association 
#'   plot for, as ENSEMBL or HGNC identifier. Character.
#' @param signal  
#' @param chrom Optional character specifying the chromosome.
#' @param pos_start Optional start position of the region for which to draw the 
#'   plot. Integer.
#' @param pos_end Optional end position of the region for which to draw the 
#'   plot. Integer.
#' @param pos Optional position for which to define the region around. Integer.
#' @param hgncid Optional HGNC identifier of the gene to define the region 
#'   around. Character.
#' @param ensemblid Optional ENSEMBL gene identifier to define the region 
#'   around. Character.
#' @param rs dbSNP RS identifier of the variant to define the region around.
#'   Character.
#' @param surround Distance around the gene to include in the region for 
#'   plotting. Default is \code{500000}.
#' @param maf_ge Numeric filtering threshold. Minor alleles with a frequency 
#'   greater than or equal to \code{maf_ge} will be included in the plot.
#' @param rsq_ge Numeric filtering threshold, specifying imputation R-squared 
#'   greater than or equal to \code{rsq_ge}.
#' @param emac_ge Numeric filtering threshold. Minor allele count must be 
#'   greater than or equal to \code{emac_ge}.
#' @param case_emac_ge Numeric filtering threshold. Case minor allele count must
#'   be greater than or equal to \code{case_emac_ge}.
#' @param priorsd Optional. Default is \code{1}. Numeric.
#' @param priorc Optional. Default is \code{1e-5}. Numeric.
#' @param cs_size Default is \code{0.95}. Numeric.
#' @param plot_ymax Numeric maximum y-axis value. Default is \code{300}.
#' @param style Character vector specifying the style or styles of plot to 
#'   output. Can be one or more of \code{c('none', 'signals', 'signal', 
#'   'classic', 'ld')}. Default is \code{'signals'}.
#' @param protein_coding_only Whether to restrict annotation to protein coding
#'   genes only. Default is \code{TRUE}.
#' @param highlight_style Highlighting style for plotting. Default is 
#'   \code{'circle'}.
#' @param dbc A database connection. Default is 
#'   \code{getOption("gtx.dbConnection", NULL)}.
#' @return Regional association plots and a \code{\linkS4class{regionplot}} 
#'   object containing the data used to draw the regional association plot.
#'
#' @details \code{regionplot()} draws a regional association plot by querying the 
#'   required data via a database connection (see \code{\link{gtxdbcheck}}). 
#'   
#'   The region to plot results over can be specified in several different ways.
#'   The region can be supplied as physical coordinates using the arguments 
#'   \code{chrom}, \code{pos}, \code{pos_start} and \code{pos_end}. 
#'   Alternatively, the region can be centered on a gene of interest, using 
#'   either the \code{hgncid} or \code{ensemblid} argument. The size of the 
#'   region around the gene can be modified using the \code{surround} argument.
#'   
#' @author Toby Johnson \email{Toby.x.Johnson@@gsk.com}
#' @family \code{\link{regionplot.data}, \link{getDataForRegionplot}}
#' @seealso \code{\link{getDataForRegionplot()}} to get the data for 
#'   \code{regionplot()}.
#' 
#' @export
regionplot <- function(analysis, 
                       entity,
                       signal, 
                       chrom, 
                       pos_start, 
                       pos_end, 
                       pos,  
                       hgncid, 
                       ensemblid, 
                       rs, 
                       surround = 500000, 
                       maf_ge, 
                       rsq_ge, 
                       emac_ge, 
                       case_emac_ge, 
                       priorsd = 1, 
                       priorc = 1e-5, 
                       cs_size = 0.95, 
                       plot_ymax = 300, 
                       style = 'signals',
                       protein_coding_only = TRUE, 
                       highlight_style = 'circle', 
                       dbc = getOption("gtx.dbConnection", NULL)) {

  if (isTRUE(getOption('gtx.debug'))) {
	futile.logger::flog.threshold(DEBUG)
  }

  # check database connection
  gtxdbcheck(dbc)

  # if optional signal argument used, check other arguments for compatibility
  if (!missing(signal)) {
    if (!(missing(maf_ge) && missing(rsq_ge) && missing(emac_ge) && missing(case_emac_ge))) {
      gtx_warn('Variant filtering will have NO EFFECT for conditional results (signal = {signal})')
    }
    if ('signals' %in% tolower(style)) {
      gtx_warn('style = \'signals\' incompatible with conditional results (signal = {signal}), using style = \'signal\' instead')
      style[tolower(style) == 'signals'] <- 'signal' 
    }
  }
  
  # Get the data for processing
  dataForRegionPlot <- getDataForRegionplot(analysis = analysis,
                                            entity = entity, 
                                            signal = signal,
                                            chrom = chrom, 
                                            pos_start = pos_start, 
                                            pos_end = pos_end, 
                                            pos = pos, 
                                            hgncid = hgncid, 
                                            ensemblid = ensemblid, 
                                            rs = rs, 
                                            surround = surround, 
                                            style = style,
                                            priorsd = priorsd, 
                                            priorc = priorc, 
                                            cs_size = cs_size, 
                                            maf_ge, 
                                            rsq_ge,
                                            emac_ge,
                                            case_emac_ge, 
                                            dbc = dbc)
  
  # Data processing for signals
  if ('signals' %in% tolower(style)) {
    dataForRegionplot@styleSignals <- styleSignals(
      pValues = dataForRegionplot@pValues,
      chrom = dataForRegionplot@xRegion$chrom,
      pos_start = dataForRegionplot@xRegion$pos_start,
      pos_end = dataForRegionplot@xRegion$pos_end,
      priorsd = priorsd,
      priorc = priorc,
      cs_size = cs_size
    )
  }
  
  # Data processing for signal
  if ('signal' %in% tolower(style)) {
    dataForRegionplot@styleSignal <- styleSignal(
      pValues = dataForRegionplot@pValues,
      priorsd = priorsd,
      priorc = priorc,
      cs_size = cs_size
      )
  }
  
  if ('ld' %in% tolower(style)) {
    dataForRegionplot@styleLd <- styleLd(
      pValues = dataForRegionplot@pValues,
      chrom = dataForRegionplot@xRegion$chrom,
      pos_start = dataForRegionplot@xRegion$pos_start,
      pos_end = dataForRegionplot@xRegion$pos_end,
      pos = pos,
      rs = rs,
      dbc = dbc
    )
  }
  
  ## Use na.rm in min, even though these should not be present,
  ##   because we want regionplot() to have robust behaviour
  ## Note, regionplot.data SQL includes pval IS NOT NULL clause
  ##   (whcih removes NULL but not NaN)
  ## Also, regionplot.data removes and warns if there are
  ##   non-finite pval
  pmin <- max(min(pvals$pval, na.rm = TRUE), 10^-plot_ymax)

  ## set par to ensure large enough margins, internal xaxs, regular yaxs,
  ## should apply for all plots generated
  oldpar <- par(mar = pmax(c(4, 4, 4, 4) + 0.1, par('mar')), 
                xaxs = 'i', yaxs = 'i') # xaxs='i' stops R expanding x axis; manual *1.04 for yaxs in regionplot.new

  if ('classic' %in% tolower(style)) {
    gl <- regionplot.new(chrom = xregion$chrom, pos_start = xregion$pos_start, pos_end = xregion$pos_end,
                 pmin = pmin, 
                 main = main, fdesc = fdesc, 
                 protein_coding_only = protein_coding_only, 
                 dbc = dbc)
    if (!missing(signal)) mtext(paste0('conditional for signal #', signal), 2, 2)
    ## best order for plotting
    ## Plot all variants with VEP annotation as blue diamonds in top layer
    pvals <- pvals[order(!is.na(pvals$impact), -log10(pvals$pval)), ]
    ## Next statement, captures the actual yvalues used when plotting
    pvals <- within(pvals, ploty <- regionplot.points(pos, pval,
                                  pch = ifelse(!is.na(impact), 23, 21),
                                  col = ifelse(!is.na(impact), rgb(0, 0, 1, .75), rgb(.33, .33, .33, .5)),
                                  bg = ifelse(!is.na(impact), rgb(.5, .5, 1, .75), rgb(.67, .67, .67, .5))))
    regionplot.highlight(gp, highlight_style = highlight_style)
  }
  if ('ld' %in% tolower(style)) {
    gl <- regionplot.new(chrom = xregion$chrom, pos_start = xregion$pos_start, pos_end = xregion$pos_end,
                 pmin = pmin, 
                 main = main, fdesc = fdesc, 
                 protein_coding_only = protein_coding_only, 
                 dbc = dbc)
    if (!missing(signal)) mtext(paste0('conditional for signal #', signal), 2, 2) # DRY should this go into regionplot.new?
    ## best order for plotting
    ## Plot variants shaded by LD
    pvals <- pvals[order(!is.na(pvals$impact), -log10(pvals$pval)), ]
    ## current db format means r<0.01 is not included in table, hence substitute missing with zero
    pvals$r[is.na(pvals$r)] <- 0. #FIXME dangerous to return incorrect numerical data
    pvals$alpha <- .5 + .25*pvals$r^2 # currently use a single alpha for col and bg for coding and noncoding, so precalculate
    with(pvals, regionplot.points(pos, pval,
                                  pch = ifelse(!is.na(impact), 23, 21),
                                  col = ifelse(!is.na(impact), rgb(0, 0, 0, alpha), rgb(.33, .33, .33, alpha)),
                                  bg = ifelse(!is.na(impact), 
                                              rgb((1 - r^2)*.67, (1 - r^2)*.67, .33*r^2 + .67, alpha),
                                              rgb(.33*r^2 + .67, (1 - r^2)*.67, (1 - r^2)*.67, alpha))))
    ## (re)draw the index variant for pairwise LD over the top, in larger size and with no transparency
    pval1 <- merge(pvals, attr(pvals, 'index_ld')[ , c('chrom', 'pos', 'ref', 'alt')])
    if (!is.null(pval1) && !is.na(pval1) && nrow(pval1) > 0) {
      with(pval1, {
        text((xregion$pos_start + xregion$pos_end)*0.5, 
              .75*gl$yline[5] + .25*gl$yline[4],
            labels = glue('Pairwise LD with index chr{chrom}:{pos}:{ref}>{alt}'),
            pos = 3, cex = 0.75)
        arrows(x0 = (xregion$pos_start + xregion$pos_end)*0.5, 
               y0 = .75*gl$yline[5] + .25*gl$yline[4], 
               y1 = .5*gl$yline[5] + .5*gl$yline[4], 
               length = 0, lty = 'dotted')
        arrows(x0 = (xregion$pos_start + xregion$pos_end)*0.5, 
               y0 = .5*gl$yline[5] + .5*gl$yline[4],
               x1 = pos, 
               y1 = .25*gl$yline[5] + .75*gl$yline[4],
               length = 0, lty = 'dotted')
        arrows(x0 = pos, 
               y0 = .25*gl$yline[5] + .75*gl$yline[4],
               y1 = pmin(-log10(pval), plot_ymax), # FIXME still not quite right
               length = 0, lty = 'dotted')
        regionplot.points(pos, pval, 
                           cex = 1.5, 
                           pch = ifelse(!is.na(impact), 23, 21),
                           col = ifelse(!is.na(impact), rgb(0, 0, 0, alpha = 1), rgb(.33, .33, .33, alpha = 1)),
                           bg = ifelse(!is.na(impact), 
                                       rgb(0, 0, 1, alpha = 1),
                                       rgb(1, 0, 0, alpha = 1)))
      })
    }

    regionplot.highlight(gp, highlight_style = highlight_style)
  }
  if ('signal' %in% tolower(style)) {
    gl <- regionplot.new(chrom = xregion$chrom, pos_start = xregion$pos_start, pos_end = xregion$pos_end,
                 pmin = pmin, 
                 main = main, fdesc = fdesc, 
                 protein_coding_only = protein_coding_only, 
                 dbc = dbc)
    if (!missing(signal)) mtext(paste0('conditional for signal #', signal), 2, 2) # DRY should this go into regionplot.new?
    ## best order for plotting is highest pp on top, with impact on top of that
    pvals <- pvals[order(!is.na(pvals$impact), pvals$pp_signal), ]
    pvals <- within(pvals, {
        pprel <- pp_signal/max(pp_signal, na.rm = TRUE) # relative posterior prob for colouring
        alpha <- .5 + .25*pprel # currently use a single alpha for col and bg for coding and noncoding, so precalculate
    })
    with(pvals, regionplot.points(pos, pval,
                                  pch = ifelse(!is.na(impact), 23, 21),
                                  cex = ifelse(cs_signal, 1.25, .75), 
                                  col = ifelse(!is.na(impact), rgb(0, 0, 0, alpha), rgb(.33, .33, .33, alpha)),
                                  bg = ifelse(cs_signal,
                                              ifelse(!is.na(impact),
                                                     rgb((1 - pprel)*.67, (1 - pprel)*.67, .33*pprel + .67, alpha),
                                                     rgb(.33*pprel + .67, (1 - pprel)*.67, (1 - pprel)*.67, alpha)),
                                              rgb(.67, .67, .67, .5))))
    regionplot.highlight(gp, highlight_style = highlight_style)
    # Add extra annotation to indicate CLEO conditional signal is being plotted
    if (!missing(signal)) {
      with(pvals[which.max(pvals$pp_signal), ], {
        #arrows(x0 = pos, y0 = .75*gl$yline[5] + .25*gl$yline[4], 
        #       y1 = -log10(pval), 
        #       length = 0, lty = 'dotted') 
        text(pos, -log10(max(pval, pmin)), labels = paste0('#', signal), pos = 3, cex = 0.75)
        # FIXME this may not work when yaxis truncation is on
      })
    }
  }
  if ('signals' %in% tolower(style)) {
    gl <- regionplot.new(chrom = xregion$chrom, pos_start = xregion$pos_start, pos_end = xregion$pos_end,
                 pmin = pmin, 
                 main = main, fdesc = fdesc, 
                 protein_coding_only = protein_coding_only, 
                 dbc = dbc)

    ## catch situation where CLEO results are present but no signals
    if ('signal' %in% names(pvals)) {
      signals <- sort(unique(na.omit(pvals$signal)))
    } else {
      signals <- NULL
    }
    if (length(signals) == 0L) {
      futile.logger::flog.warn('No CLEO results, using single signal results instead')
      ## no signals, so use single signal results as if they were the CLEO results
      pvals <- within(pvals, {
        cs_cleo <- cs_signal
        pp_cleo <- pp_signal
        signal <- ifelse(cs_signal, 'default', NA)
      })
      signals <- 'default'
    }
    
    ## best order for plotting
    ## Plot variants coloured/sized by credible set, with VEP annotation is diamond shape in top layer
    pvals <- pvals[order(!is.na(pvals$impact), pvals$pp_cleo), ]

    ## colvec is vector of identifying colour for each signal
    if (length(signals) == 1L) {
        colvec <- rgb(1, 0, 0) # red
    } else {
        colvec <- rainbow(length(signals), start = .6667, end = .3333) # blue-magenta-red-yellow-green avoiding cyan part of spectrum
    }

    bg_cleo <- rep(rgb(.67, .67, .67, alpha = .5), nrow(pvals))
    for (idx in 1:length(signals)) {
        ww <- which(pvals$signal == signals[idx])
        if (length(ww) > 0L) {
            pprel <- pvals$pp_cleo[ww]
            pprel <- pprel/max(pprel, na.rm = TRUE) # relative pp *within signal*
            bg_cleo[ww] <- rgb((col2rgb(colvec[idx])[1]/255 - .67)*pprel + .67, 
                               (col2rgb(colvec[idx])[2]/255 - .67)*pprel + .67,
                               (col2rgb(colvec[idx])[3]/255 - .67)*pprel + .67, 
                               alpha = .5)
        }
    }

    with(pvals, regionplot.points(pos, pval,
                                  pch = ifelse(!is.na(impact), 23, 21), 
                                  cex = ifelse(cs_cleo, 1.25, .75),
                                  bg = bg_cleo, 
                                  col = ifelse(!is.na(impact), rgb(0, 0, 0, .5), rgb(.33, .33, .33, .5))))

    ## legend indicating signals
    ## ssi = summary stats for index variants
    ssi <- attr(pvals, 'index_cleo')
    if (!is.null(ssi) && nrow(ssi) > 1) {
      ssi <- ssi[match(signals, ssi$signal), ] # order to match signals and colvec
      ssi$keypos <- seq(from = xregion$pos_start, to = xregion$pos_end, length.out = nrow(ssi) + 2)[rank(ssi$pos) + 1]
      arrows(x0 = ssi$keypos, y0 = .75*gl$yline[5] + .25*gl$yline[4], 
             y1 = .5*gl$yline[5] + .5*gl$yline[4], 
             length = 0, lty = 'dotted')
      arrows(x0 = ssi$keypos, y0 = .5*gl$yline[5] + .5*gl$yline[4],
             x1 = ssi$pos, y1 = .25*gl$yline[5] + .75*gl$yline[4],
             length = 0, lty = 'dotted')
      arrows(x0 = ssi$pos, y0 = .25*gl$yline[5] + .75*gl$yline[4],
             y1 = pmin(-log10(ssi$pval), plot_ymax), # FIXME still not quite right
             length = 0, lty = 'dotted')
      points(ssi$keypos, rep(.75*gl$yline[5] + .25*gl$yline[4], nrow(ssi)),
             pch = 21, col = rgb(.33, .33, .33, .5), bg = colvec, cex = 1)
      text(ssi$keypos, rep(.75*gl$yline[5] + .25*gl$yline[4], nrow(ssi)),
           labels = paste0('#', signals), pos = 4, cex = 0.75)
      text(mean(ssi$keypos), .75*gl$yline[5] + .25*gl$yline[4],
           labels = 'CLEO index variants', pos = 3, cex = 0.75)
      #legend("bottomleft", pch = 21, ,
      #       pt.bg = colvec, legend=paste0('#', signals),
      #       horiz=T, bty="n", cex=.5)
    } else {
      text((xregion$pos_start + xregion$pos_end)/2., .75*gl$yline[5] + .25*gl$yline[4],
           labels = 'No CLEO index variants', pos = 3, cex = 0.75)
    }
    regionplot.highlight(gp, highlight_style = highlight_style)
  } 

  ## restore par
  par(oldpar)

  ## for when called from within GUI tools, return 
  ## the pvalue dataframe.  Should also return information on the genelayout. 
  return(invisible(pvals))
}

#' Get data for a regional association plot.
#' 
#' \code{getDataForRegionplot()} returns all of the data required to draw a
#'   regional association plot with the \code{regionplot()} function.
#' 
#' @inheritParams regionplot
#' @return A \code{\linkS4class{regionplot}} object containing the data needed 
#'   to draw a regional association plot.
#'
#' @family \code{\link{regionplot}} functions.
#' @seealso \code{\link{regionplot()}} to get the data and draw regional 
#'   association plots in one step.
#' 
#' @export
getDataForRegionplot <- function(analysis, 
                                 entity, 
                                 signal, 
                                 chrom, 
                                 pos_start, 
                                 pos_end, 
                                 pos, 
                                 hgncid, 
                                 ensemblid, 
                                 rs, 
                                 surround = 500000,
                                 style = 'signals',
                                 priorsd = 1,
                                 priorc = 1e-5,
                                 cs_size = 0.95, 
                                 maf_ge,
                                 rsq_ge,
                                 emac_ge, 
                                 case_emac_ge, 
                                 dbc = getOption("gtx.dbConnection", NULL)) {
  
  # Check the database connection
  gtxdbcheck(dbc)
  
  # Initialise a new regionplot object
  dataForRegionplot <- new('regionplot')
  
  # Determine x-axis range from arguments
  dataForRegionplot@xRegion <- gtxregion(chrom = chrom, 
                                         pos_start = pos_start, 
                                         pos_end = pos_end, 
                                         pos = pos, 
                                         hgncid = hgncid, 
                                         ensemblid = ensemblid, 
                                         rs = rs, 
                                         signal = signal, 
                                         analysis = analysis, 
                                         entity = entity, 
                                         surround = surround,
                                         dbc = dbc)
  
  # If required, determine entity and associated info including entity_label
  dataForRegionplot@xEntity <- gtxentity(analysis = analysis, 
                                         entity = entity, 
                                         hgncid = hgncid, 
                                         ensemblid = ensemblid)
  
  # Get the p-values
  dataForRegionplot@pValues <- getPValues(analysis = analysis,
                                          signal = signal,
                                          chrom = dataForRegionplot@xRegion$chrom,
                                          pos_start = dataForRegionplot@xRegion$pos_start,
                                          pos_end = dataForRegionplot@xRegion$pos_end,
                                          style = style,
                                          maf_ge = maf_ge,
                                          rsq_ge = rsq_ge,
                                          emac_ge = emac_ge,
                                          case_emac_ge = case_emac_ge,
                                          entityWhere = dataForRegionplot@xEntity$entity_where,
                                          xRegionLabel = dataForRegionplot@xRegion$label,
                                          dbc = dbc)
  
  # Set the title of the plot (used to be 'main')
  dataForRegionplot@plotTitle <- gtxanalysis_label(analysis = analysis, 
                                                   entity = dataForRegionplot@xEntity, 
                                                   signal = signal, 
                                                   nlabel = TRUE, 
                                                   dbc = dbc)
  
  # Set the subtitle of the plot (used to be 'fdesc')
  # in future we may need to pass maf_lt and rsq_lt as well  
  if (missing(signal)) {
    dataForRegionplot@plotSubTitle <- gtxfilter_label(maf_ge = maf_ge, 
                                                      rsq_ge = rsq_ge,
                                                      emac_ge = emac_ge, 
                                                      case_emac_ge = case_emac_ge, 
                                                      analysis = analysis)
  } else {
    dataForRegionplot@plotSubTitle <- 'All variants (after pre-CLEO filtering)'
  }
  
  # Get index SNPs to highlight if the pos or rs arg was used
  if (!missing(pos)) {
    dataForRegionplot@indexSNPs <- dataForRegionplot@pValues[dataForRegionplot@pValues$pos %in% pos, c('pos', 'pval')]
  }
  if (!missing(rs)) {
    
    # Query this rs id
    rsQuery <- getDataFromDB(connectionType = 'SQL',
                             connectionArguments = list(dbc,
                                                        sprintf('SELECT chrom, pos, ref, alt 
                                                                FROM sites 
                                                                WHERE %s;',
                                                                gtxwhere(chrom = dataForRegionplot@xRegion$chrom, 
                                                                         rs = rs)),
                                                        uniq = FALSE, 
                                                        zrok = TRUE))
    dataForRegionplot@indexSNPs <- rbind(dataForRegionplot@indexSNPs, 
                                         merge(dataForRegionplot@pValues, 
                                               rsQuery, 
                                               all.x = FALSE, 
                                               all.y = TRUE)[ , c('pos', 'pval')])
  }
  
  return(dataForRegionplot)
}

# getPValues gets the p-values for regionplot
# The arguments for this function are almost all inherited from 
# getDataForRegionplot
# entityWhere is the entity for plotting, which can be produced with gtxentity()
# xRegionLabel is the label for the region to plot, which can be produced with
# gtxregion()
getPValues <- function(analysis,
                       signal, 
                       chrom,
                       pos_start,
                       pos_end,
                       style = 'signals',
                       maf_ge, 
                       rsq_ge,
                       emac_ge, 
                       case_emac_ge, 
                       entityWhere,
                       xRegionLabel,
                       dbc = getOption("gtx.dbConnection", NULL)) {
  
  # Query marginal p-values
  if (missing(signal)) {
    
    gtx_info('Querying marginal p-values')
    t0 <- as.double(Sys.time())
    
    pValues <- getDataFromDB(connectionType = 'SQL',
                             connectionArguments = list(dbc, 
                                                        sprintf('SELECT gwas_results.chrom, gwas_results.pos, gwas_results.ref, gwas_results.alt, pval, impact %s 
                                                                FROM %sgwas_results 
                                                                LEFT JOIN vep 
                                                                USING (chrom, pos, ref, alt) 
                                                                WHERE %s AND %s AND %s AND %s AND pval IS NOT NULL;',
                                                                if (any(c('signal', 'signals') %in% tolower(style))) ', beta, se, rsq, freq' else '', 
                                                                gtxanalysisdb(analysis), 
                                                                gtxwhat(analysis1 = analysis),
                                                                entityWhere, # (entity=...) or (True)
                                                                gtxwhere(chrom, 
                                                                         pos_ge = pos_start, 
                                                                         pos_le = pos_end, 
                                                                         tablename = 'gwas_results'),
                                                                gtxfilter(maf_ge = maf_ge, 
                                                                          rsq_ge = rsq_ge, 
                                                                          emac_ge = emac_ge, 
                                                                          case_emac_ge = case_emac_ge, 
                                                                          analysis = analysis)),
                                                        uniq = FALSE))
    
    t1 <- as.double(Sys.time())
    gtx_info('Query returned {nrow(pValues)} variants in query region {xRegionLabel} in {round(t1 - t0, 3)}s.')
    
  } else {
    
    gtx_info('Querying conditional p-values')
    t0 <- as.double(Sys.time())
    
    pValues <- getDataFromDB(connectionType = 'SQL', 
                             connectionArguments = list(dbc, 
                                                        sprintf('SELECT gwas_results_cond.chrom, gwas_results_cond.pos, gwas_results_cond.ref, gwas_results_cond.alt, 
                                                                beta_cond AS beta, se_cond AS se, impact 
                                                                FROM %sgwas_results_cond 
                                                                LEFT JOIN vep 
                                                                USING (chrom, pos, ref, alt) 
                                                                WHERE %s AND %s AND %s;',
                                                                gtxanalysisdb(analysis), 
                                                                where_from(analysisu = analysis, 
                                                                           signalu = signal, 
                                                                           tablename = 'gwas_results_cond'), 
                                                                entityWhere, # (entity=...) or (True)
                                                                gtxwhere(chrom, 
                                                                         pos_ge = pos_start, 
                                                                         pos_le = pos_end, 
                                                                         tablename = 'gwas_results_cond')), 
                                                        uniq = FALSE))
    
    t1 <- as.double(Sys.time())
    gtx_info('Query returned {nrow(pvals)} variants in query region {xRegionLabel} in {round(t1 - t0, 3)}s.')
    
    pValues$pval <- with(pValues, pchisq((beta/se)^2, df = 1, lower.tail = FALSE))
  }
  
  if (any(!is.finite(pValues$pval))) {
    gtx_warn('Removing {sum(!is.finite(pValues$pval))} variants with non-finite P-values')
    
    pValues <- pValues[is.finite(pValues$pval), ]
    
    gtx_debug('After non-finite P-values removed, {nrow(pValues)} variants in query region {xRegionLabel}')
  }
  
  pValues <- within(pValues, impact[impact == ''] <- NA)
  pValues <- pValues[order(pValues$pval), ]
  
  return(pValues)
}

# Data processing for the 'signals' style
# Returns a named list
styleSignals <- function(pValues,
                         analysis,
                         chrom,
                         pos_start,
                         pos_end,
                         priorsd,
                         priorc,
                         cs_size){
  
  ## The code below is from regionplot.data()
  
  futile.logger::flog.debug('Finemapping under single signal assumption')
  
  signalsPValues <- fm_signal(pValues, 
                              priorsd = priorsd, 
                              priorc = priorc, 
                              cs_size = cs_size,
                              cs_only = FALSE)
  
  # get CLEO credible sets only, since we effectively left join with this, 
  # for plot colouring
  fmResults <- fm_cleo(analysis = analysis, 
                       chrom = chrom, 
                       pos_start = pos_start, 
                       pos_end = pos_end,
                       priorsd = priorsd, 
                       priorc = priorc, 
                       cs_size = cs_size, 
                       cs_only = TRUE)
    
  # note fm_cleo already prints logging messages about number of signals
    
  if (nrow(fmResults) > 0) {
    
    # sort by decreasing posterior probability
    # match each row or pvals with *first* match in fmResults, thus linking any 
    # variants in more than one credible set, with the one for the signal for 
    # which it has higher probability
    
    fmResults <- fmResults[order(fmResults$pp_cleo, decreasing = TRUE), ]
    
    signalsPValues <- cbind(signalsPValues,
                            fmResults[match(with(signalsPValues, 
                                                 paste(chrom, pos, ref, alt,
                                                       sep = '_')), 
                                            with(fmResults, 
                                                 paste(chrom, pos, ref, alt, 
                                                       sep = '_'))),
                                      c('signal', 'cs_cleo', 'pp_cleo')])
    
    ## FIXME in case one variant is in more than one credible
    ## set, it would be better to cbind like we do above with signal,
    ## but then cbind with the *marginal* pp's from aggregating over signals
    
    # Convert NA values of pp_cleo to zero to make them safe with sorting and 
    # plotting
    signalsPValues$pp_cleo[is.na(signalsPValues$pp_cleo)] <- 0. 
    
    # Set NA values of cs_cleo to FALSE to be safe with tables and plotting
    signalsPValues$cs_cleo[is.na(signalsPValues$cs_cleo)] <- FALSE 
  }
  
  styleSignalsObj <- new('styleSignals')
  styleSignalsObj@signalsPValues <- signalsPValues
  styleSignalsObj@index_cleo <- attr(fmResults, "index_cleo")
  
  return(styleSignalsObj)
}

# Data processing for the 'signal' style
# Returns a data frame with the p-values for the 'signal' style
styleSignal <- function(pValues, priorsd, priorc, cs_size){
  
  # The code below is from regionplot.data()
  
  futile.logger::flog.debug('Finemapping under single signal assumption')
  signalPValues <- fm_signal(pValues, 
                             priorsd = priorsd, 
                             priorc = priorc, 
                             cs_size = cs_size, 
                             cs_only = FALSE)
  
  # Remove unused attributes from df
  attr(signalPValues, "params_signal") <- NULL
  attr(signalPValues, "nullpp_signal") <- NULL
  
  styleSignalObj <- new('styleSignal')
  styleSignalObj@signalPValues <- signalPValues
  
  return(styleSignalObj)
}

# Data processing for the 'classic' style
styleClassic <- function(){
  
}

# Data processing for the 'ld' style
styleLd <- function(pValues, 
                    chrom, 
                    pos_start, 
                    pos_end, 
                    pos, 
                    rs, 
                    dbc = getOption('gtx.dbConnection', NULL)){
  
  ## The code below is from regionplot.data()
  
  
  # Merge with LD information in userspace code since we need to pull the 
  # p-values first to find the top hit or otherwise chosen variant
  
  # Create a variable to store index variant chrom/pos/ref/alt
  indexVariant <- NULL
  
  # If a single pos argument was provided, try to use this as index variant
  if (!missing(pos)) {
    
    if (identical(length(pos), 1L)) {
      indexVariant <- pValues[pValues$pos == pos, c('chrom', 'pos', 'ref', 'alt'), drop = FALSE]
      
      if (identical(nrow(indexVariant), 1L)) {
        indexVariant$r <- 1
        gtx_debug('Selected pairwise LD index chr{indexVariant$chrom}:{indexVariant$pos}:{indexVariant$ref}>{indexVariant$alt} 
                  (selected by pos argument)')
      } else {
        gtx_warn('Skipping pos argument [ {paste(pos, collapse = \', \')} ] for 
                 pairwise LD index selection because {nrow(indexVariant)} 
                 variants match')
        indexVariant <- NULL
      }
    } else {
      gtx_warn('Skipping pos argument [ {paste(pos, collapse = \', \')} ] for 
               pairwise LD index selection because multiple values')
    }
  }
  
  # If a single rs argument was provided, try to use this as index variant
  if (is.null(indexVariant) && !missing(rs)) {
    
    if (identical(length(rs), 1L)) {
      rsQuery <- getDataFromDB(connectionType = 'SQL',
                               connectionArguments = list(dbc,
                                                          sprintf('SELECT chrom, pos, ref, alt 
                                                                  FROM sites 
                                                                  WHERE %s;',
                                                                  gtxwhere(chrom = chrom, 
                                                                           rs = rs)),
                                                          uniq = FALSE, 
                                                          zrok = TRUE))
      
      indexVariant <- merge(pValues, rsQuery)[ , c('chrom', 'pos', 'ref', 'alt'), 
                                               drop = FALSE]
      
      if (identical(nrow(indexVariant), 1L)) {
        indexVariant$r <- 1
        gtx_debug('Selected pairwise LD index 
                  chr{indexVariant$chrom}:{indexVariant$pos}:{indexVariant$ref}>{indexVariant$alt} 
                  (selected by rs argument)')
        
      } else {
        gtx_warn('Skipping rs argument [ {paste(rs, collapse = \', \')} ] for 
                 pairwise LD index selection because {nrow(indexVariant)} 
                 variants match')
        indexVariant <- NULL
      }
    } else {
      gtx_warn('Skipping rs argument [ {paste(rs, collapse = \', \')} ] for 
               pairwise LD index selection because multiple values')            
    }
  }
  
  # Next, pick variant with at least one ld value and with smallest P-value 
  if (is.null(indexVariant)) {
    
    # Query all variants present on LHS (chrom1, pos1, ref1, alt1) of pairwise 
    # LD table
    # Note that if db guarantees 1:1 between variants in the pairwise TABLE ld, 
    # and the per-variant TABLE ldref, this can be done more efficiently
    # Notes:  cannot use gtxwhere because of nonstandard names chrom1, pos1
    hasLd <- getDataFromDB(connectionType = 'SQL',
                           connectionArguments = list(dbc, 
                                                      sprintf('SELECT chrom1 AS chrom, pos1 AS pos, ref1 AS ref, alt1 AS alt, True AS has_ld \
                                                              FROM ld \
                                                              WHERE chrom1 = \'%s\' AND pos1 >= %s AND pos1 <= %s \
                                                              GROUP BY chrom1, pos1, ref1, alt1;',
                                                              sanitize1(chrom, 
                                                                        values = c(1:22, 'X')),
                                                              sanitize1(pos_start, 
                                                                        type = 'int'),
                                                              sanitize1(pos_end, 
                                                                        type = 'int')),
                                                      uniq = FALSE, 
                                                      zrok = TRUE))
    
    hasLd <- within(merge(hasLd, pValues[ , c('chrom', 'pos', 'ref', 'alt', 'pval'), drop = FALSE], 
                          all.x = FALSE, 
                          all.y = TRUE), 
                    has_ld[is.na(has_ld)] <- FALSE)
    
    if (nrow(hasLd) == 0L || all(!hasLd$has_ld)) {
      gtx_warn('Skipping pairwise LD index selection because no variants have 
               pairwise LD data')
    } else {
      
      # order so we can warn how many smaller P-value variants were skipped 
      hasLd <- hasLd[order(hasLd$pval), , drop = FALSE]
      
      if (hasLd$has_ld[1]) {
        # smallest P-value variant has LD, so okay to use this
        indexVariant <- hasLd[1, c('chrom', 'pos', 'ref', 'alt'), drop = FALSE]
        indexVariant$r <- 1
        gtx_debug('Selected pairwise LD index chr{pval1$chrom}:{pval1$pos}:{pval1$ref}>{pval1$alt} 
                  (smallest P-value)')
      } else {
        whichHasLd <- which(hasLd$has_ld)[1]
        indexVariant <- hasLd[whichHasLd, c('chrom', 'pos', 'ref', 'alt'), drop = FALSE]
        indexVariant$r <- 1
        gtx_warn('Pairwise LD index selection skipped {whichHasLd - 1} variants 
                 with -log10(p)<={round(log10(hasLd$pval[whichHasLd]/hasLd$pval[1]), 2)} 
                 smaller, with no pairwise LD data')
        gtx_debug('Selected pairwise LD index chr{indexVariant$chrom}:{indexVariant$pos}:{indexVariant$ref}>{indexVariant$alt} 
                  (smallest P-value with pairwise LD data)')
      }
    }
  }
  
  if (!is.null(indexVariant)) {
    
    # Stop if hasLd was all FALSE or had zero rows
    stopifnot(identical(nrow(indexVariant), 1L))
    
    gtx_debug('Querying pairwise LD with index chr{indexVariant$chrom}:{indexVariant$pos}:{indexVariant$ref}>{indexVariant$alt}')
    
    pairwiseLd <- getDataFromDB(connectionType = 'SQL',
                         connectionArguments = list(dbc, 
                                                    sprintf('SELECT chrom2 AS chrom, pos2 AS pos, ref2 AS ref, alt2 AS alt, r FROM ld
                                                            WHERE chrom1=\'%s\' AND pos1=%s AND ref1=\'%s\' AND alt1=\'%s\';',
                                                            indexVariant$chrom, 
                                                            indexVariant$pos, 
                                                            indexVariant$ref, 
                                                            indexVariant$alt), # should we sanitize
                                                    uniq = FALSE, 
                                                    zrok = TRUE))
    
    ldPValues <- merge(pValues, rbind(indexVariant, pairwiseLd), all.x = TRUE, all.y = FALSE)
    
    # sort by decreasing r^2 (will be resorted for plotting, so makes 
    # a .data() call work as a useful proxy search
    ldPValues <- ldPValues[order(ldPValues$r^2, decreasing = TRUE), ]
    
    # threshold seems high but pairwise LD often missing form low frequency 
    # variants
    if (mean(is.na(ldPValues$r)) > 0.25) {
      gtx_warn('Pairwise LD with index 
               chr{pval1$chrom}:{pval1$pos}:{pval1$ref}>{pval1$alt} missing for 
               {round(mean(is.na(pvals$r))*100)}% of variants, check warnings 
               and debugging messages for possible cause')
    }
    
  } else {
    gtx_warn('Selection of pairwise LD index failed, check warnings and 
             debugging messages for possible cause')
    ldPValues$r <- NA
  }
  
  styleLdObj <- new('styleLd')
  styleLdObj@ldPvalues <- ldPValues
  styleLdObj@indexVariant <- indexVariant
  
  return(styleLdObj)
}

#' @export
regionplot.data <- function(analysis, entity, signal, 
                            chrom, pos_start, pos_end, pos, 
                            hgncid, ensemblid, rs, 
                            surround = 500000,
                            style = 'signals',
                            priorsd = 1, priorc = 1e-5, cs_size = 0.95, 
                            maf_ge, rsq_ge,
                            emac_ge, case_emac_ge, 
                            dbc = getOption("gtx.dbConnection", NULL)) {
    ## check database connection
    gtxdbcheck(dbc)

    style <- tolower(style)
    
    ## Determine x-axis range from arguments
    xregion <- gtxregion(chrom = chrom, pos_start = pos_start, pos_end = pos_end, pos = pos, 
                         hgncid = hgncid, ensemblid = ensemblid, rs = rs, 
                         signal = signal, analysis = analysis, entity = entity, 
                         surround = surround,
                         dbc = dbc)
    chrom = xregion$chrom # note, overwriting command line arguments
    pos_start = xregion$pos_start
    pos_end = xregion$pos_end
    
    ## If required, determine entity and associated info including entity_label
    xentity <- gtxentity(analysis, entity = entity, hgncid = hgncid, ensemblid = ensemblid)
    
    ## always query marginal p-values
    ## seems more flexible to query CLEO results separately and merge within R code
    if (missing(signal)) {
      gtx_info('Querying marginal p-values')
      t0 <- as.double(Sys.time())
      pvals <- getDataFromDB(connectionType = 'SQL',
                             connectionArguments = list(dbc, 
                                                        sprintf('SELECT gwas_results.chrom, gwas_results.pos, gwas_results.ref, gwas_results.alt, pval, impact %s 
                                                                FROM %sgwas_results 
                                                                LEFT JOIN vep 
                                                                USING (chrom, pos, ref, alt) 
                                                                WHERE %s AND %s AND %s AND %s AND pval IS NOT NULL;',
                                                                if (any(c('signal', 'signals') %in% tolower(style))) ', beta, se, rsq, freq' else '', 
                                                                gtxanalysisdb(analysis), 
                                                                gtxwhat(analysis1 = analysis),
                                                                xentity$entity_where, # (entity=...) or (True)
                                                                gtxwhere(chrom, 
                                                                         pos_ge = pos_start, 
                                                                         pos_le = pos_end, 
                                                                         tablename = 'gwas_results'),
                                                                gtxfilter(maf_ge = maf_ge, 
                                                                          rsq_ge = rsq_ge, 
                                                                          emac_ge = emac_ge, 
                                                                          case_emac_ge = case_emac_ge, 
                                                                          analysis = analysis)),
                                                        uniq = FALSE))
      t1 <- as.double(Sys.time())
      gtx_info('Query returned {nrow(pvals)} variants in query region {xregion$label} in {round(t1 - t0, 3)}s.')
    } else {
      gtx_info('Querying conditional p-values')
      t0 <- as.double(Sys.time())
      pvals <- getDataFromDB(connectionType = 'SQL', 
                             connectionArguments = list(dbc, 
                                                        sprintf('SELECT gwas_results_cond.chrom, gwas_results_cond.pos, gwas_results_cond.ref, gwas_results_cond.alt, 
                                                                beta_cond AS beta, se_cond AS se, impact 
                                                                FROM %sgwas_results_cond 
                                                                LEFT JOIN vep 
                                                                USING (chrom, pos, ref, alt) 
                                                                WHERE %s AND %s AND %s;',
                                                                gtxanalysisdb(analysis), 
                                                                where_from(analysisu = analysis, 
                                                                           signalu = signal, 
                                                                           tablename = 'gwas_results_cond'), 
                                                                xentity$entity_where, # (entity=...) or (True)
                                                                gtxwhere(chrom, 
                                                                         pos_ge = pos_start, 
                                                                         pos_le = pos_end, 
                                                                         tablename = 'gwas_results_cond')), 
                                                        uniq = FALSE))
      t1 <- as.double(Sys.time())
      gtx_info('Query returned {nrow(pvals)} variants in query region {xregion$label} in {round(t1 - t0, 3)}s.')
      pvals$pval <- with(pvals, pchisq((beta/se)^2, df = 1, lower.tail = FALSE))
    }
      
    if (any(!is.finite(pvals$pval))) {
      gtx_warn('Removing {sum(!is.finite(pvals$pval))} variants with non-finite P-values')
      pvals <- pvals[is.finite(pvals$pval), ]      
      gtx_debug('After non-finite P-values removed, {nrow(pvals)} variants in query region {xregion$label}')
    }
    
    if (any(c('signal', 'signals') %in% tolower(style))) {
        futile.logger::flog.debug('Finemapping under single signal assumption')
        ## cs_only = FALSE since we still want to plot/return variants not in the credible set
        pvals <- fm_signal(pvals, priorsd = priorsd, priorc = priorc, cs_size = cs_size, cs_only = FALSE)
    }
    
    if ('signals' %in% tolower(style)) {
        ## get CLEO credible sets only, since we effectively left join with this, for plot colouring
        fmres <- fm_cleo(analysis = analysis, chrom = chrom, pos_start = pos_start, pos_end = pos_end,
                         priorsd = priorsd, priorc = priorc, cs_size = cs_size, cs_only = TRUE)
        ## note fm_cleo already prints logging messages about number of signals
        if (nrow(fmres) > 0) {
            ## sort by decreasing posterior probability, match each
            ## row or pvals with *first* match in fmres
            ## thus linking any variants in more than one
            ## credible set, with the one for the signal for which
            ## it has higher probability
            fmres <- fmres[order(fmres$pp_cleo, decreasing = TRUE), ]
            pvals <- cbind(pvals,
                           fmres[match(with(pvals, paste(chrom, pos, ref, alt, sep = '_')),
                                       with(fmres, paste(chrom, pos, ref, alt, sep = '_'))),
                                 c('signal', 'cs_cleo', 'pp_cleo')])
            ## FIXME in case one variant is in more than one credible
            ## set, it would be better to cbind like we do above with signal,
            ## but then cbind with the *marginal* pp's from aggregating over signals
            pvals$pp_cleo[is.na(pvals$pp_cleo)] <- 0. # make zero to be safe with sorting/plotting
            pvals$cs_cleo[is.na(pvals$cs_cleo)] <- FALSE # make FALSE to be safe with tables/plotting
            ## Propagate CLEO index variants as attr()ibute
            attr(pvals, 'index_cleo') <- attr(fmres, 'index_cleo')
        } else {
            # do nothing
        }
    }

    if ('ld' %in% tolower(style)) {
        # merge with LD information in userspace code since we need to
        # pull the p-values first to find the top hit or otherwise chosen variant

        # note if style='ld', sort order is by ld with index variant,
        # otherwise sort by pval (see else block)

        # Determine Index Variant
        pval1 <- NULL # will be used to store chrom/pos/ref/alt of the index variant
                      # and using NULL to indicate not successfully selected
        # First, if a single pos argument was provided, try to use this as index variant
        if (!missing(pos)) {
          if (identical(length(pos), 1L)) {
            ## funny syntax to avoid subset(pvals, pos %in% pos)
            pval1 <- pvals[pvals$pos == pos, c('chrom', 'pos', 'ref', 'alt'), drop = FALSE]
            if (identical(nrow(pval1), 1L)) {
              pval1$r <- 1
              gtx_debug('Selected pairwise LD index chr{pval1$chrom}:{pval1$pos}:{pval1$ref}>{pval1$alt} (selected by pos argument)')
            } else {
              gtx_warn('Skipping pos argument [ {paste(pos, collapse = \', \')} ] for pairwise LD index selection because {nrow(pval1)} variants match')
              pval1 <- NULL
            }
          } else {
            gtx_warn('Skipping pos argument [ {paste(pos, collapse = \', \')} ] for pairwise LD index selection because multiple values')
          }
        }
        # Next, if a single rs argument was provided, try to use this as index variant
        if (is.null(pval1) && !missing(rs)) {
          if (identical(length(rs), 1L)) {
            # query this rs id
            qrs <- getDataFromDB(connectionType = 'SQL',
                                 connectionArguments = list(dbc,
                                                            sprintf('SELECT chrom, pos, ref, alt 
                                                                FROM sites 
                                                                WHERE %s;',
                                                                    gtxwhere(chrom = chrom, 
                                                                             rs = rs)),
                                                            uniq = FALSE, 
                                                            zrok = TRUE))
            # use merge as a way to subset on match by all of chrom/pos/ref/alt
            # note default merge is all.x = FALSE, all.y = FALSE
            pval1 <- merge(pvals, qrs)[ , c('chrom', 'pos', 'ref', 'alt'), drop = FALSE]
            if (identical(nrow(pval1), 1L)) {
              pval1$r <- 1
              gtx_debug('Selected pairwise LD index chr{pval1$chrom}:{pval1$pos}:{pval1$ref}>{pval1$alt} (selected by rs argument)')
            } else {
              gtx_warn('Skipping rs argument [ {paste(rs, collapse = \', \')} ] for pairwise LD index selection because {nrow(pval1)} variants match')
              pval1 <- NULL
            }
          } else {
            gtx_warn('Skipping rs argument [ {paste(rs, collapse = \', \')} ] for pairwise LD index selection because multiple values')            
          }
        }
        # Next, pick variant with at least one ld value and with smallest P-value 
        if (is.null(pval1)) {
          # Query all variants present on LHS (chrom1, pos1, ref1, alt1) of pairwise LD table
          # Note that if db guarantees 1:1 between variants in the pairwise TABLE ld, and the
          # per-variant TABLE ldref, this can be done more efficiently
          # Notes:  cannot use gtxwhere because of nonstandard names chrom1, pos1
          has_ld <- getDataFromDB(connectionType = 'SQL',
                                  connectionArguments = list(dbc, 
                                                             sprintf('SELECT chrom1 AS chrom, pos1 AS pos, ref1 AS ref, alt1 AS alt, True AS has_ld \
                                                                 FROM ld \
                                                                 WHERE chrom1 = \'%s\' AND pos1 >= %s AND pos1 <= %s \
                                                                 GROUP BY chrom1, pos1, ref1, alt1;',
                                                                     sanitize1(xregion$chrom, 
                                                                               values = c(1:22, 'X')), # should be a type for chrom  
                                                                     sanitize1(xregion$pos_start, 
                                                                               type = 'int'),
                                                                     sanitize1(xregion$pos_end, 
                                                                               type = 'int')),
                                                             uniq = FALSE, 
                                                             zrok = TRUE))
          has_ld <- within(merge(has_ld, pvals[ , c('chrom', 'pos', 'ref', 'alt', 'pval'), drop = FALSE],
                                 all.x = FALSE, all.y = TRUE),
                           has_ld[is.na(has_ld)] <- FALSE)
          if (nrow(has_ld) == 0L || all(!has_ld$has_ld)) {
            gtx_warn('Skipping pairwise LD index selection because no variants have pairwise LD data')
          } else {
            # order so we can warn how many smaller P-value variants were skippedcount 
            has_ld <- has_ld[order(has_ld$pval), , drop = FALSE]
            if (has_ld$has_ld[1]) {
              # smallest P-value variant has LD, so okay to use this
              pval1 <- has_ld[1, c('chrom', 'pos', 'ref', 'alt'), drop = FALSE]
              pval1$r <- 1
              gtx_debug('Selected pairwise LD index chr{pval1$chrom}:{pval1$pos}:{pval1$ref}>{pval1$alt} (smallest P-value)')
            } else {
              w_has_ld <- which(has_ld$has_ld)[1] # guaranteed length >=1 by checks above
              pval1 <- has_ld[w_has_ld, c('chrom', 'pos', 'ref', 'alt'), drop = FALSE]
              pval1$r <- 1
              gtx_warn('Pairwise LD index selection skipped {w_has_ld - 1} variants with -log10(p)<={round(log10(has_ld$pval[w_has_ld]/has_ld$pval[1]), 2)} smaller, with no pairwise LD data')
              gtx_debug('Selected pairwise LD index chr{pval1$chrom}:{pval1$pos}:{pval1$ref}>{pval1$alt} (smallest P-value with pairwise LD data)')
            }
          }
        }
        if (!is.null(pval1)) {
          stopifnot(identical(nrow(pval1), 1L)) # This will fail if has_ld was all FALSE or had zero rows, see above
          gtx_debug('Querying pairwise LD with index chr{pval1$chrom}:{pval1$pos}:{pval1$ref}>{pval1$alt}')
          ld1 <- getDataFromDB(connectionType = 'SQL',
                               connectionArguments = list(dbc, 
                                                          sprintf('SELECT chrom2 AS chrom, pos2 AS pos, ref2 AS ref, alt2 AS alt, r 
                                                                    FROM ld
                                                                    WHERE chrom1=\'%s\' AND pos1=%s AND ref1=\'%s\' AND alt1=\'%s\';',
                                                                  pval1$chrom, 
                                                                  pval1$pos, 
                                                                  pval1$ref, 
                                                                  pval1$alt), # should we sanitize
                                                          uniq = FALSE, 
                                                          zrok = TRUE))
          pvals <- merge(pvals, rbind(pval1, ld1), all.x = TRUE, all.y = FALSE)
          # sort by decreasing r^2 (will be resorted for plotting, so makes 
          # a .data() call work as a useful proxy search
          pvals <- pvals[order(pvals$r^2, decreasing = TRUE), ]
          attr(pvals, 'index_ld') <- pval1
          if (mean(is.na(pvals$r)) > 0.25) { # threshold seems high but pairwise LD often missing form low frequency variants
            gtx_warn('Pairwise LD with index chr{pval1$chrom}:{pval1$pos}:{pval1$ref}>{pval1$alt} missing for {round(mean(is.na(pvals$r))*100)}% of variants, check warnings and debugging messages for possible cause')
          }
        } else {
          gtx_warn('Selection of pairwise LD index failed, check warnings and debugging messages for possible cause')
          pvals$r <- NA
          attr(pvals, 'index_ld') <- NA
        }
    } else {
        ## sort by increasing pval
        pvals <- pvals[order(pvals$pval), ]
    }
    
    attr(pvals, 'analysis') <- analysis
    attr(pvals, 'region') <- xregion
    attr(pvals, 'entity') <- xentity
    return(pvals)
}

#' @export
regionplot.new <- function(chrom, pos_start, pos_end, pos, 
                           hgncid, ensemblid, rs, surround = 500000, 
                           pmin = 1e-10, main, fdesc, 
                           protein_coding_only = TRUE,   
                           dbc = getOption("gtx.dbConnection", NULL)) {
  gtxdbcheck(dbc)

  ## Determine x-axis range from arguments
  xregion <- gtxregion(chrom = chrom, pos_start = pos_start, pos_end = pos_end, pos = pos, 
                       hgncid = hgncid, ensemblid = ensemblid, rs = rs, surround = surround,
                       dbc = dbc)
  chrom = xregion$chrom
  pos_start = xregion$pos_start
  pos_end = xregion$pos_end

  ## Determine y-axis upper limit
  ## removed +0.5 since space above allocated by regionplot.genelayout()
  ymax <- ceiling(-log10(pmin))

  ## Determine amount of y-axis space needed for gene annotation
  gl <- regionplot.genelayout(chrom, pos_start, pos_end, ymax, protein_coding_only = protein_coding_only)
  gtx_debug('gene layout set up for ymax={gl$ymax} yline[1]={gl$yline[1]} yline[4]={gl$yline[4]} yline[5]={gl$yline[5]}')
  ## Set up plotting area
  plot.new()
  plot.window(c(pos_start, pos_end), range(gl$yline * 1.04, na.rm=TRUE)) # y axis *1.04 adjustment for using par(yaxs='i')

  ## Draw axes, and axis labels
  abline(h = 0, col = "grey")
  with(list(xpretty = pretty(c(pos_start, pos_end)*1e-6)),
       axis(1, at = xpretty*1e6, labels = xpretty))
  with(list(ypretty = pretty(c(0, ymax))),
       axis(2, at = ypretty[ypretty <= ymax], las = 1))
  mtext(paste0("chr", chrom, " position (Mb)"), 1, 3)
  mtext(expression(paste("Association ", -log[10](paste(P, "-value")))), 2, 3)

  ## Add recombination rate and gene annotation
  ##  regionplot.recombination determines position range from par("usr")
  regionplot.recombination(chrom, yoff = mean(gl$yline[2:3]))
  ##  regionplot.genedraw uses previously determined layout
  regionplot.genedraw(gl)

  ## Add title, scaled to fit if necessary, CHANGE TO USE mtext.fit() FIXME
  if (!missing(main)) {
    xplt <- par("plt")[2] - par("plt")[1] # figure as fraction of plot, assumes no subsequent changes to par("mar")
    mtext(main, 3, 1,
          cex = min(1., xplt/strwidth(main, units = 'figure')))
  }
  if (!missing(fdesc)) {
      mtext(fdesc, 3, 0, cex = 0.5)
  }

  ## Draw box last to overdraw any edge marks
  box()

  return(invisible(gl))
}

#' @export
regionplot.genedraw <- function(gl) {
  arrows(gl$genelayout$pos_start, gl$genelayout$y, x1 = gl$genelayout$pos_end,
         length = 0, lwd = 2, col = "blue")
  if (length(gl$genelayout$label) > 0) {
    text(gl$genelayout$pos_mid, gl$genelayout$y - .25*strheight("M", cex = gl$cex),
         gl$genelayout$label,
         adj = c(.5, 1), font = 3, cex = gl$cex, col = "blue")
  }
  return(invisible(NULL))
}

# each gene is assigned to a line (iline) such that the horizontal line
# and label will not overlap neighboring genes.  The y positions are
# then evenly distributed [add more explanation]
#
# returns a list with elements:
# cex:  value of cex used when computing the layout
# ymax: value of ymax used when computing the layout
# yline:  vector of length 5 delimiting the regions for the genes track, 
#                   zero [min for points], max for points, max for whole plot   
#' @export
regionplot.genelayout <- function (chrom, pos_start, pos_end, ymax, cex = 0.75, 
				   protein_coding_only = TRUE, 
                                   dbc = getOption("gtx.dbConnection", NULL)) {
  gtxdbcheck(dbc)
  
  yplt <- par("plt")[4] - par("plt")[3] # figure as fraction of plot, assumes no subsequent changes to par("mar")
  xplt <- par("plt")[2] - par("plt")[1]
  xusr <- pos_end - pos_start # par("usr")[2] - par("usr")[1] 
  return(with(getDataFromDB(connectionType = 'SQL',
                            connectionArguments = list(getOption('gtx.dbConnection_cache_genes', dbc),
                                                       sprintf('SELECT min(pos_start) AS pos_start, max(pos_end) AS pos_end, hgncid, ensemblid 
                                                               FROM genes 
                                                               WHERE %s %s 
                                                               GROUP BY ensemblid, hgncid 
                                                               ORDER BY pos_start', 
                                                               gtxwhere(chrom = chrom, 
                                                                        pos_end_ge = pos_start, 
                                                                        pos_start_le = pos_end),
                                                               if (protein_coding_only) 'AND genetype=\'protein_coding\'' else ''),
                                                       uniq = FALSE, 
                                                       zrok = TRUE)),
              {
                label <- ifelse(hgncid != '', as.character(hgncid), as.character(ensemblid)) # could also check NULL or NA?
                ## compute start and end plot positions for each gene, using larger of transcript line and gene name annotation
                pos_mid <- .5*(pos_start + pos_end)
                xwidth <- strwidth(label, units = 'figure', cex = cex)*xusr/xplt
                xpad <- max(10000, strwidth('M', units = 'figure', cex = cex)*xusr/xplt)
                xstart <- pmin(pos_start, pos_mid - .5*xwidth) - xpad
                xend <- pmax(pos_end, pos_mid + .5*xwidth) + xpad
                ## layout engine based on putting each gene on the highest line possible
                iline <- rep(NA, length(label))
                if (length(label) > 0) {
                  iline[1] <- 0
                  if (length(label) > 1) {
                    for (idx in 2:length(label)) {
                      iline[idx] <- with(rbind(aggregate(xend[1:(idx - 1)], by = list(tline = iline[1:(idx - 1)]), FUN = max),
                                               data.frame(tline = max(iline[1:(idx - 1)]) + 1, x = -Inf)),
                                         min(tline[x < xstart[idx]]))
                    }
                  }
                }
                fy <- 0 # fraction of total figure y space needed
                if (length(iline) > 0) {
                  fy <- max(iline + 1)*2*strheight("M", units = "figure", cex = cex)/yplt # fraction of total figure region needed
                  if (fy > 0.5) {
                    futile.logger::flog.warn('Squashing gene annotation to fit within .5 of plot area')
                    fy <- 0.5
                  }
                }
                yd1 = ymax/(0.8-fy) * 0.1
                yd2 = ymax/(0.8-fy) * fy
                yline <- c(-(yd1 + yd2), -yd1, 0, ymax, ymax + yd1)
                y <- numeric(0)
                if (length(iline) > 0) {
                  y <- -(iline/max(iline + 1)*yd2 + yd1)
                }
                list(cex = cex, ymax = ymax, yline = yline, 
                     genelayout = data.frame(label, pos_start, pos_end, pos_mid, iline, y))
              }))
}

regionplot.points <- function(pos, pval, labels, 
                              pch = 21, bg = rgb(.67, .67, .67, .5), col = rgb(.33, .33, .33, .5), cex = 1, 
                              pos.labels = 3, col.labels = rgb(0, 0, 0), cex.labels = 0.5, 
                              ymax, 
                              suppressWarning = FALSE) {
  # auto-detecting ymax (==yline[4]), works because:
  # par('usr')[5] == yline[5] # true since changed to par(yaxs='i')
  # par('usr')[4] == yline[1] # ibid
  # yline[5]-yline[4] is 10% of plot always (by genelayout())
  # reverse the *1.04 applied in regionplot.new plus
  # plus a 0.5 'epsilon' to avoid effect of numerical imprecision on floor()  
  if (missing(ymax)) ymax <- floor(sum(par('usr')[c(3, 4)]*c(.1, .9)/1.04) + 0.5)
  y <- -log10(pval)
  f <- y > ymax
  if (missing(labels)) {
    points(pos, ifelse(f, ymax, y), pch = pch, col = col, bg = bg, cex = cex)
  } else {
    text(x = pos, y = ifelse(f, ymax, y), labels = labels, 
         pos = pos.labels, col = col.labels, cex = cex.labels)
  }
  if (any(f) && !suppressWarning) {
    # would be nice to more cleanly overwrite y axis label
    axis(2, at = ymax, labels = substitute({}>=ymax, list(ymax = ymax)), las = 1)
    # draw at left edge, otherwise multiple regionplot.points() generate 
    # partly overlapping  (and hence illegible) warnings
    text(par('usr')[1] + strwidth('M', cex = 0.5), 
         ymax, '(plot truncated)', adj = c(0, 0.5), cex = 0.5)
    return(FALSE)
  }
  return(ifelse(f, ymax, y))
}

regionplot.highlight <- function(pvals, highlight_style) {
  stopifnot(is.data.frame(pvals))
  if (all(c('pos', 'pval') %in% names(pvals))) {
        pvals <- subset(pvals, !is.na(pos) & !is.na(pval))
        if (nrow(pvals) > 0) {
            if ('circle' %in% tolower(highlight_style)) {
                regionplot.points(pvals$pos, pvals$pval, 
                                  pch = 1, cex = 2, col = rgb(0, 0, 0, .75),
                                  suppressWarning = TRUE) # suppress since points already drawn
            }
            if ('pos' %in% tolower(highlight_style)) {
                regionplot.points(pvals$pos, pvals$pval, 
                                  labels = pvals$pos,
                                  suppressWarning = TRUE) # suppress since points already drawn
            }
            ## in future will add options if 'rs' %in% highlight_style etc.
        }
    } else {
        # silently do nothing, e.g. if pvals==data.frame(NULL)
    }
    return(invisible(NULL))
}

regionplot.recombination <- function(chrom, pos_start, pos_end, yoff = -.5, 
                                     dbc = getOption("gtx.dbConnection", NULL)) {
  gtxdbcheck(dbc)

  ## Recombination segments (r) that fall wholly or partly within [pos_start, pos_end]
  ## have r.pos_end >= pos_start and r.pos_start <= pos_end
  if (missing(pos_start)) pos_start = floor(par("usr")[1])
  if (missing(pos_end)) pos_end = ceiling(par("usr")[2])

  with(getDataFromDB(connectionType = 'SQL',
                     connectionArguments = list(dbc, 
                                                sprintf('SELECT pos_start, pos_end, recombination_rate 
                                                        FROM genetic_map 
                                                        WHERE %s', 
                                                        gtxwhere(chrom = chrom, 
                                                                 pos_end_ge = pos_start,
                                                                 pos_start_le = pos_end)),
                                                uniq = FALSE, 
                                                zrok = TRUE)),
       {
         abline(h = yoff, col = "grey")
         yscale <- (par("usr")[4] - yoff)*.75/max(recombination_rate)
         lines(c(pos_start, pos_end[length(pos_end)]),
               yoff + c(recombination_rate, recombination_rate[length(recombination_rate)])*yscale,
               type = "s", col = "cyan3")
         with(list(ypretty = pretty(c(0, max(recombination_rate)))),
              axis(4, at = yoff + ypretty*yscale, labels = ypretty, las = 1))
       })
  mtext("Recombination rate (cM/Mb)", 4, 3)
  return(invisible(NULL))
}

#library for region plot
#Author: Li Li, modified based on Toby's code



#draw recomination rate 
#' @export
recomb.draw <- function(chr, from.bp, to.bp, ylo, yhi, recomb){
  if (is.integer(chr)) chr <- as.character(chr)
  if(chr == "23") chr <- "X"
  chr <- ifelse(substr(chr, 1, 3) == "chr", chr, paste("chr", chr, sep = ""))
  #recomb <- read.table(paste("genetic_map_chr", chr, ".txt", sep=""), header=T)
  #load recombination info
  recomb.chr <- recomb[[chr]][-1]
  idx.start<- which(recomb.chr[,1] < from.bp)
  idx.end<- which(recomb.chr[,1] > to.bp)
  keep.recomb <-recomb.chr[idx.start[length(idx.start)]: idx.end[1],]
  big.range <- yhi -ylo
  max.curr<- max(keep.recomb[,2], na.rm = T)
  #print (paste("Max recomb rate in the region is", max.curr))
  ticks.curr <- pretty(c(0, max.curr))
  max.curr<- ticks.curr[length(ticks.curr )]
  lines(keep.recomb[,1], ylo + ( ( keep.recomb[,2] / max.curr ) * big.range), type="l", 
        col="lightblue", lwd=1)
  axis(4, at=ylo + big.range /max.curr * ticks.curr, labels=ticks.curr, las=1, hadj = 0)
  mtext("Recomb rate (cM/Mb)", side=4, at=(ylo+big.range/2), line=2.5) 
}



################################################################################
#gwas: gwas results with columns for MARKER, chr, POS (bp), qcflag, genoflag, pvalues
#chr pos: hit snp
#flanking: flanking rang (kb) for plotting
#col.pvals: column names for p values to be plotted by different shapes (21:25, 3,4,7:14) 
#plabel: label for p value column
#col.r2: column name for r2 to index SNP. Background colored by LD to index SNP (white to red)
################################################################################
#' @export
regionplot.multiP <- function(gwas, chr, pos, flanking = 250, col.r2=NULL,
                              col.pvals = c("P1","P2"), plabel = NULL, gencode = gencode, recomb=recomb,
                              main = "", main.sub = NULL, miny = 8, crity = -log10(5e-08)) {
  to.bp <- min(pos+flanking *1000, max(gwas$POS, na.rm=T))
  from.bp <- max(pos - flanking *1000, min(gwas$POS, na.rm = T))
  stopifnot(to.bp > from.bp)
  stopifnot(all(col.pvals %in% names(gwas)))
  sel <- gwas$CHROM == chr & gwas$POS >= from.bp & gwas$POS <= to.bp
  cat(sum(sel), "SNPs selected\n")
  pos.sel<- gwas$POS[sel]
  minP <- 1e-50 #minimum value for p value
  for( c in col.pvals) gwas[[c]][!is.na(gwas[[c]]) & gwas[[c]] < minP] <- minP
  lp <- -log10(as.matrix(subset(gwas, sel, col.pvals))) # lp = -log10(p)
  if (is.null(crity)) crity <- -log10(0.05/nrow(lp))
  if (is.null(plabel)) plabel <- col.pvals
  r2.sel <- rep(1, sum(sel))
  if(!is.null(col.r2)) {
    r2.sel<- gwas[[col.r2]][sel]
    r2.sel[is.na(r2.sel)]<- 0
    r2.sel <- pmin(1, pmax(0, 1-r2.sel))
  }
 
 
  plot.new()
  yhi<- max(c(miny, 1.05 * max(lp, na.rm = TRUE)))
  y.tick<- pretty(c(0, yhi))
  yhi <- y.tick[length(y.tick)]
  scale  <- yhi/8
  plot.window(c(from.bp, to.bp), c(-6 *scale,  yhi))
  
  abline(h = seq(0, 7, 1), col = "lightgrey", lty = "dotted")
  abline(h = crity, col = "red", lty = "dashed")
  ppch <- c(21:25, 3,4,7:14) #rep(21:25, 5)
  index <- pos.sel == pos 
  for (idx in 1:ncol(lp)) {
    points(pos.sel[-index], lp[-index , idx], pch = ppch[idx], col = "black",  bg = rgb(1, r2.sel[-index], r2.sel[-index]), cex = 1)
    if(any(index)) #index snp
      points(pos, lp[index, idx],  pch = ppch[idx], col = rgb(0, 0, 1),  bg = rgb(0, 0, 1), cex = 1.3)
  }   
  legend("topright",   ncol = 2,   bty = "n",  pch = ppch[1:ncol(lp)], 
           text.col =  "black", cex = 0.8,  legend = plabel)
  for(i in 0:10){
    polygon(from.bp + (c(0, 1, 1, 0)+i) *(to.bp-from.bp)/80, yhi- c(0.5, 0.5, 1, 1)*scale, 
            border = "grey", col = rgb(1, 1-i/10, 1-i/10)) 
  }
  text(from.bp, yhi - 0 *scale , expression(paste("LD Map Type:","r"^"2")), adj = 0, cex = 0.8)
  text(from.bp+seq(0.5, 10.5, 2) * (to.bp-from.bp)/80, rep(yhi, 6)- 1.5*scale, 
       labels =c(0, 0.2, 0.4, 0.6, 0.8,1),  cex = 0.7)
  
  xm <- pretty(c(from.bp, to.bp)*1e-6) # x marks
  axis(1, at = xm*1e6, labels = xm)
  axis(2, at = y.tick, las = 1)
  mtext(expression(-log[10](italic(P))), side=2, at=yhi/2, line=2.5) 
  recomb.draw(chr, from.bp, to.bp, -1*scale, yhi, recomb)                        
  gene.draw(paste("chr", chr, sep = ""), from.bp, to.bp, gencode, 
            yhi = -2*scale, ylo = -6*scale,exony = 0.05*scale, genecex = 0.7)
  mtext(paste("Chromosome ", chr, " genomic position (Mb)", sep = ""), side = 1, line = 2.5)
  title(main = main)  # ylab = expression(-log[10](italic(P))),
  if(!is.null(main.sub))
    title(main = main.sub, cex.main = 1, col.main = "blue", line = 0.5)
  box()
  return(invisible(NULL))
}


