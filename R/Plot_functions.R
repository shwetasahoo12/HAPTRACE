#### Plot functions ####

#' Plot for the QTL effect
#'
#' @description
#' Generates the density plot for the QTL effects.
#'
#'
#' @param QTL_eff_file QTL file with information on QTL effects.
#' @param color color of the density plot.
#'
#' @returns shows density plot depicting the distribution of the QTL effects.
#' @seealso [HAPTRACE::QTLeffects()]
#' @export
#'
#' @importFrom stats density
#' @importFrom graphics polygon
#'
#' @examples example_QTL <- QTLeffects(n_QTL_Chr = 1, nChr = 5, len_Chr = c(10, 8, 6, 4, 2),
#' shape_gamma = 0.4)
#' plot <- QTL_densityplot(QTL_eff_file = example_QTL, color = "forestgreen")
QTL_densityplot = function(QTL_eff_file, color = "steelblue") {
  if (missing(QTL_eff_file)) {
    stop("QTL effect file is not defined")
  }
  else {
    QTLefffile <- stats::density(as.numeric(QTL_eff_file[, 'Effect']))
    Plot <- plot(QTLefffile,
                 frame = FALSE,
                 col = color,
                 main = "Density plot of QTL effect")
    graphics::polygon(QTLefffile, col = color)
  }
}


#' Plot for True Breed Composition
#'
#' @description
#' Generates barplot of true breed composition of the animals with haplotype information
#' coded according to the breed of origin.
#'
#' @param trueBC True breed composition estimated by the haplotype information coded
#' according to the breed of origin.
#' @param color Colors specifying the number of populations contributing to the true breed composition.
#' @param legendlabels Population label to add as a legend to the plot.
#'
#' @returns Returns an admixture plot of animals where x axis shows number of animals
#' and y axis shows the proportion of different breeds in different colors.
#' @seealso [HAPTRACE::trueBreedComp()]
#' @export
#'
#' @importFrom graphics barplot
#'
#'
#'
#' @examples popAdm <- matrix(c(-1, 2, 3, -4, 2,
#' 1,-1,2,2,3,
#' 1,2,-3,-4,-5,
#' 2,-3,2,4,-5), byrow = TRUE, nrow = 4, ncol = 5)
#' TBC <- trueBreedComp(Haplotype = popAdm)
#' Admixplot(trueBC = TBC, color = c("#FFC107", "#2E7D32", "#FB8072", "#386CB0", "#4d4a49"),
#' legendlabels = c("A", "B", "C", "D", "E"))
Admixplot = function(trueBC, color, legendlabels)
{
  if (missing(trueBC)) {
    stop(
      "Define the true breed composition estimated from coded haplotype. Refer trueBreedComp function."
    )
  }
  else if (missing(color)) {
    stop("Please define the number of animals to be simulated per generation")
  }
  else if (missing(legendlabels)) {
    stop(
      "Please define the labels of the population present while estimating true breed composition"
    )
  }
  else if (length(legendlabels) != nrow(trueBC)) {
    stop("Check the number of populations and define the labels")
  }
  else if (length(color) != nrow(trueBC)) {
    stop("The number of colors should match the number of populations")
  }
  else {
    xLabels <- c(1:ncol(trueBC))
    graphics::barplot(
      trueBC,
      col = color,
      xlab = "Animals",
      ylab = "Ancestry",
      border = NA,
      names.arg = xLabels
    )
    # Add the legend
    graphics::legend(
      "top",
      inset = c(-0.15, -0.15),
      legend = legendlabels,
      pch = 15,
      col = color,
      horiz = TRUE,
      xpd = TRUE,
      bty = "n"
    )
  }
}


#' Plot for recombination points and recombination probability
#'
#' @description
#' Function to generate plot for recombination points and recombination probability.
#'
#' @param recList List of recombination points with each list containing information
#' of recombination points of an animal.
#' @param col Color of the recombination probability plot.
#' @param map Map file of the markers for the population.
#' @param n_markers Total number of markers.
#'
#' @returns This function returns the plots for recombination points and recombination probability.
#' @seealso [HAPTRACE::generateMAP()], [HAPTRACE::RecombinationPoint()]
#'
#' @export
#'
#' @importFrom graphics axis abline plot clip par barplot
#'
#' @examples map_rec <- generateMAP(nChr = 2, n_markers = 5, len_Chr = 10)
#' poprec <- matrix(data = sample(c(0,1), 1000, replace = TRUE), nrow = 100, ncol = 10)
#' RecPoint <- RecombinationPoint(map = map_rec, Haplotype = poprec, recL = 2, nChr = 2)
#' reclist <- RecPoint$recList
#' recom_plot <- recProbplot(recList = reclist, n_markers = 10, map = map_rec)
recProbplot = function(recList, n_markers, col = "steelblue", map) {
  if (missing(recList) || !is.list(recList)) {
    stop("Please define the recombination points in list")
  }
  else if (missing(n_markers)) {
    stop("Define the number of markers")
  } else if (missing(map)) {
    stop("Define map file")
  } else{
    chr_endpoint <- as.vector(cumsum(table(map$Chr)))
    chr_swap_rm <- chr_endpoint + 1
    ##### Unlist the recombination
    list <- unlist(recList)
    zero_vector <- vector(length = n_markers)
    zero_vector[] <- 0
    # Create a table of counts
    counts <- table(list)

    # Get numeric names (the elements) and their counts
    indices <- as.numeric(names(counts))
    values <- as.vector(counts)

    # Update zero vector
    zero_vector[indices] <- values

    ## Calculate probability
    nChr <- length(as.matrix(as.numeric(unique(map$Chr))))
    prob_count <- zero_vector
    prob_count[chr_swap_rm] <- 0
    prob = prob_count / length(prob_count)
    prob_plot = graphics::plot(
      prob,
      type = "l",
      col = col,
      xlab = "" ,
      ylab = "Probability",
      main = "Recombination Probability"
    )
    chr_names = paste0("Chr ", 1:nChr)
    graphics::axis(
      1,
      at = chr_endpoint,
      labels = chr_names,
      col = "black",
      col.axis = "black",
      hadj = 1.7,
      padj = -2.4,
      las = 3,
      cex.axis = 0.9
    )
    graphics::clip(
      x1 = par("usr")[1],
      x2 = par("usr")[2],
      y1 = par("usr")[3],
      y2 = par("usr")[4]
    )
    graphics::abline(v = chr_endpoint, col = "black")

    ### recombination points
    rec_vector <- zero_vector
    rec_vector[chr_swap_rm] <- 0
    ### plot distribution of recombination points
    nchr_recptcount <- graphics::barplot(table(rec_vector), xlab = "Distribution of Recombination points")
    return(list(
      recProb = prob,
      recProbplot = prob_plot,
      recpt_count = nchr_recptcount
    ))
  }
}


#' Plot for Haplotype segments
#'
#' @description
#' Function to generate plot of haplotype segments inherited from different breeds.It
#' utilises the coded haplotype information according to the breed of origin and generates plot
#' at all chromosomes.
#'
#' @param Haplotype Haplotype information of the animals to be plotted.
#' @param colors Colors to be defined for different populations.
#' @param map Mapfile of the population comprising of marker ID, chromosome number and marker position.
#' @param nChr Haploid number of chromosomes(n).
#' @param legendlabels Population label to add as a legend to the plot.
#'
#' @returns Haplotype plot showing blocks from different population for animals.
#' @seealso [HAPTRACE::generateMAP()]
#'
#' @export
#'
#' @importFrom graphics image axis mtext legend abline
#'
#' @examples Haplotype <- matrix(data = sample(c(-1, 1, -2, 2, -3, 3), 100, replace = TRUE),
#' nrow = 10, ncol = 10)
#' Haplotype_ID <- PopulationID(Population = Haplotype, PopCode = "H1")
#' map <- generateMAP(nChr = 2, n_markers = 5, len_Chr = 5)
#' colors = c("gold", "forestgreen", "steelblue")
#' plot <- plotHaplo(Haplotype = Haplotype_ID, colors = colors, map = map,
#' nChr = 2, legendlabels = c("A", "B", "C"))
plotHaplo = function(Haplotype, colors, map, nChr, legendlabels) {
  if (missing(Haplotype)) {
    stop("Define the coded haplotype matrix")
  }
  if (is.null(rownames(Haplotype))) {
    stop("Row names should be defined with AnimalID for haplotype matrix")
  }
  if (!is.matrix(Haplotype)) {
    stop("Haplotype information is not in matrix form")
  }
  if (base::nrow(Haplotype) < 2L){
    stop("Haplotype matrix must have atleast two rows")
  }
  if (base::nrow(Haplotype) %% 2 != 0){
    stop("Haplotype matrix must have even number of rows",
         "(two per animal)")
  }
  else if (missing(colors)) {
    stop("Define the colors per population")
  }
  else if (length(colors) != length(unique(c(abs(Haplotype))))) {
    stop(paste("The number of colors (", length(colors), ") should match the number of populations (", length(unique(c(abs(Haplotype)))), ")",
               sep = ""))
  }
  else if (missing(map)) {
    stop("Define the map file of the population")
  }
  else if (missing(nChr)) {
    stop("Define the number of haploid number of chromosomes (N)")
  }
  else if (nChr > 40) {
    warning(
      "The haploid number of chromosomes exceeds 40. Please check the haploid number of chromosomes and proceed with caution!"
    )
  }
  else if (missing(legendlabels)) {
    stop("Define the labels of the population")
  }
  else if (length(legendlabels) != length(unique(c(abs(Haplotype))))) {
    stop("The number of colors should match the number of populations")
  }
  else {
    nmarkers_chr <- as.vector(cumsum(table(map$Chr)))
    x <- abs(Haplotype)
    xLabels <- c(1:ncol(x)) ## number of markers
    yLabels <- rownames(Haplotype) ## number of animals
    reverse <- nrow(x):1
    yLabels <- yLabels[reverse]
    x <- x[reverse, ]
    min <- min(x)
    max <- max(x)
    graphics::image(
      1:length(xLabels),
      1:length(yLabels),
      t(x),
      col = colors,
      xlab = "",
      ylab = "",
      axes = FALSE,
      zlim = c(min, max)
    )
    graphics::axis(
      BELOW <- 1,
      at = 1:length(xLabels),
      labels = xLabels,
      cex.axis = 0.7
    )
    graphics::axis(
      LEFT <- 2,
      at = 1:length(yLabels),
      labels = yLabels,
      las = HORIZONTAL <- 1,
      cex.axis = 0.65
    )
    # Add custom x-axis with labels
    chr_names <- paste0("Chr ", 1:nChr)
    graphics::axis(
      1,
      at = nmarkers_chr,
      labels = chr_names,
      col = "black",
      col.axis = "black",
      hadj = 1.7,
      padj = -0.5,
      las = 3,
      cex.axis = 0.9
    )
    ## Add vertical line
    graphics::abline(v = nmarkers_chr + 0.5, col = "black")
    graphics::mtext(side = 4, paste("Number of animals with haplotype: ", nrow(x) /
                                      2))
    graphics::mtext(side = 3, paste("Number of markers:", ncol(x)))
    # Add the legend
    graphics::legend(
      "top",
      inset = c(-0.15, -0.15),
      legend = legendlabels,
      pch = 15,
      col = colors,
      horiz = TRUE,
      xpd = TRUE,
      bty = "n"
    )
  }
}
