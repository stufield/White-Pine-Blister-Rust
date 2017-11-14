

#' Calculate Limited Cumulative Sum
#' 
#' Calculate the cumulative sum of a set of numbers within a vector. The
#' difference between this and \code{\link{cumsum}} is that it is a sliding window
#' approach, so sums are not necessarily calculated over the entire length of
#' the vector. When cut == length(x) then it is the same as \code{\link{cumsum}}. Also, for
#' entries < cut, i.e. the beginning, the entries returned will be identical to
#' \code{\link{cumsum}}.
#' 
#' This function is primarily used in making figures for the pine12 manuscript.
#' 
#' @param x The vector to be summed across
#' @param cut The size of the window to sum within (the cutoff)
#' @return A vector of the sums of the sliding window for the cumulative sums
#' @section Warning:
#' If the cut >= length(x), a warning is generated and \code{cumsum2} 
#' reverts to \code{\link{cumsum}}
#' @author Stu Field
#' @seealso \code{\link{cumsum}}
#' @examples
#' 
#' cumsum2(1:20, 5)
#' cumsum2(1:20, cut=20)
#' cumsum(1:20)
#' r.vec <- sample(1:20, 100, replace=TRUE) # random vector
#' cumsum2(r.vec, 5)
#' 
#' @export cumsum2
cumsum2 <- function(x, cut) {

   hi <- (cut+1):length(x)
   lo <- hi - (cut - 1)

   if ( cut >= length(x) ) {
      warning(sprintf("The cut (%i) >= length of x (%i)! Defaulting to cumsum(x)\n",
                      cut, length(x)))
      return(cumsum(x))
   } else {
      one <- cumsum(x[1:cut])
      two <- unlist(lapply(1:length(hi), function(i) sum(x[lo[i]:hi[i]])))
      #for (i in 1:length(hi)) out[cut+i] <- sum(x[lo[i]:hi[i]]) # loop option
      c(one, two)
   }
}


#' Calculate Dead Trees (pine12)
#' 
#' Calculate cumulative number of dead trees during a population projection
#' 
#' Used internally during \code{\link{pine12}}.
#' 
#' @param PM The (post) matrix containing the stages (cols) with 
#' solutions at all time steps (rows)
#' @param s Survivorship parameters (diagonal)
#' @param t Transition parameters (sub-diagonal)
#' @param cc Cost parameters (as placed in the cost matrix along the diagonal)
#' @return A matrix of the cumulative number of dead trees, rows = time steps,
#' cols = classes.
#' @author Stu Field
#' @seealso \code{\link{pine12}}
#' @references Thanks to Jesse Drendel for simplifying the code considerably.
#' @examples
#' 
#' # see pine12 output
#' 
#' @export calcDeadTrees
calcDeadTrees <- function(PM, s, t, cc) {

   D  <- matrix(0, ncol=ncol(PM), nrow=nrow(PM))
   Sv <- rep(s, 2)
   Tv <- rep(t, 2)
   for ( i in 1:nrow(PM) ) {
      for ( j in 2:ncol(PM) ) {
         D[i, j] <- PM[i, j] * (1 - cc[j] * (Sv[j] + Tv[j]))
      }
   }
   D <- as.data.frame(round(D,1))
   names(D) <- paste0(" C", 1:12)
   return(D)
}



#' Lambda Growth (pine12)
#' 
#' Determines the proportional growth (>1) or decline (<1) of the whole
#' population in a stage structured model. This is different from true
#' \eqn{lambda} in matrix modeling as it is not a property of the projection
#' matrix itself, but is calculated post-process based on the solutions.
#' 
#' @param x A numeric vector of the total population over time.
#' @return Proportional growth (x_t+1 / x_t)
#' @note Used in pine12
#' @author Stu Field
#' @seealso \code{\link{do.call}}
#' @examples
#' 
#' pop <- (1:10)^2
#' LambdaGrow(pop)
#' 
#' @export LambdaGrow
LambdaGrow <- function(x) {
   do.call(c, sapply(seq_along(x), function(i) x[i] / x[i - 1]))
}



#' Mating Subroutine (pollen)
#' 
#' Same as mating36() assuming there is a constant pollen cloud determining
#' Hardy-Weinberg ratios.
#' 
#' A sub-routine used to determine the number of offspring produced by the
#' current population according to Hardy-Weinberg predictions. In this
#' function, the mating population is the top TWO classes, Sub-adults and
#' Mature adults. The mating matrix (Y) does not explicitly have to have ONLY
#' the mating population, because other classes will be multiplied by zero. In
#' the function, the fitness and fecundity matrices are immediately multiplied
#' together entry-wise (i.e. the Hadamard product).
#' 
#' @param Y A matrix containing the mating population with the classes as rows
#' & genotypes as cols. Columns should be ordered AA, Aa, aa. Only the final
#' TWO classes (5 & 6) mate.
#' @param W The fitness matrix. A matrix of identical dimensions to Y with entries 
#' containing the combination of the fitness matrix (W) which contains the 
#' relative fitnesses of the mating population and the fecundity matrix (F) which 
#' contains the relative fecundities (offspring produced per individual per class).
#' @param FM The fecundity matrix. A matrix of identical dimensions to Y with the 
#' entries containing the relative fecundities of the mating population. 
#' Only the mating population should have non-zero entries.
#' @param a The allele frequency (i.e. q) of the recessive allele to be kept
#' constant according to the assumption of a pollen cloud determining
#' fertilization. 0 < a < 1.
#' @return A matrix of identical dimensions to Y containing non-zero entries in
#' the first row only representing the newly produced seeds in the population.
#' This matrix can then simply be added to the current, or intermediate,
#' solution to obtain the current solution.
#' @note Used in pine36() assuming a constant pollen cloud. Jesse's Function for
#' calculating Mating freqs & subsequent offspring 6 Classes; Top 2 mate 
#' Both infected & susceptible reproduce For use in pollen cloud simulation 
#' with constant freq.R
#' @author Stu Field, Jesse Drendel, Simon Tavener
#' @seealso \code{\link{mating36}}
#' @references Hardy-Weinberg
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' # See pine36()
#' 
#' @export mating36p
mating36p <- function(Y, W, FM, a) {

   WF     <- W * FM
   F_seed <- matrix(0, nrow=nrow(Y), ncol=ncol(Y))
   nR1R1  <- sum(Y[c(5, 6, 11, 12), 1])
   nR1R2  <- sum(Y[c(5, 6, 11, 12), 2])
   nR2R2  <- sum(Y[c(5, 6, 11, 12), 3])
   Total  <- nR1R1 + nR1R2 + nR2R2    # total number of reproductive adults

   if ( Total > 0 ) {
     p <- 1 - a
     q <- a
     #***************************
     fR1R1 <- sum((Y * WF)[,1])
     fR1R2 <- sum((Y * WF)[,2])
     fR2R2 <- sum((Y * WF)[,3])
     #***************************
     F.R1R1 <- p * (fR1R2/2 + fR1R1)
     F.R1R2 <- q * (fR1R2/2 + fR1R1) + p * (fR1R2/2 + fR2R2)
     F.R2R2 <- q * (fR1R2/2 + fR2R2)
     #***************************
     F_seed[1, ] <- c(F.R1R1, F.R1R2, F.R2R2)
   }

   return(F_seed)
}



#' Mating Frequency Subroutine (pine36)
#' 
#' A sub-routine used to determine the number of offspring produced by the
#' current population according to Hardy-Weinberg predictions, for the
#' pine36 function.
#' 
#' In this function, the mating population is the top TWO classes, sub-adults
#' and mature adults. The mating matrix (Y) does not explicitly have to have
#' ONLY the mating population, because other classes will be multiplied by
#' zero. In the function, the fitness and fecundity matrices are immediately
#' multiplied together entry-wise (i.e. the Hadamard product).
#' 
#' @param Y The mating populaton matrix. A matrix containing the mating 
#' population with the classes as rows and genotypes as cols. Columns should be 
#' ordered AA, Aa, aa. Only the final \emph{2} classes (5 & 6) reproduce.
#' @param W The fitness matrix. A matrix of identical dimensions to Y with entries 
#' containing the combination of the fitness matrix (W) which contains the 
#' \emph{relative} fitnesses of the mating population and the fecundity matrix (F) 
#' which contains the relative fecundities (offspring produced per individual per class).
#' @param FM The fecundity matrix. A matrix of identical dimensions to Y with the 
#' entries containing the relative fecundities of the mating population. 
#' Only the mating population should have non-zero entries.
#' @return A matrix of identical dimensions to Y containing non-zero entries in
#' the first row only representing the newly produced seeds in the population.
#' This matrix can then simply be added to the current, or intermediate,
#' solution to obtain the current solution.
#' @note Used in pine36
#' @author Stu Field, Jesse Drendel, Mike Antolin, Simon Tavener
#' @seealso \code{\link{pine36}}, \code{\link{mating36}}
#' @references Jesse Drendel
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' # See pine36()
#' 
#' @export mating36
mating36 <- function(Y, W, FM) {

   F_seed <- matrix(0, nrow=nrow(Y), ncol=ncol(Y))
   WF     <- W * FM    # Hadamard product
   nR1R1  <- sum(Y[c(5, 6, 11, 12), 1])
   nR1R2  <- sum(Y[c(5, 6, 11, 12), 2])
   nR2R2  <- sum(Y[c(5, 6, 11, 12), 3])
   Total  <- nR1R1 + nR1R2 + nR2R2    # total number of reproductive adults

   if ( Total > 0 ) {
     p      <- (nR1R2 + 2 * nR1R1) / (2 * Total)
     q      <- (nR1R2 + 2 * nR2R2) / (2 * Total)
     fR1R1  <- sum((Y * WF)[, 1])
     fR1R2  <- sum((Y * WF)[, 2])
     fR2R2  <- sum((Y * WF)[, 3])
     F.R1R1 <- p * (fR1R2/2 + fR1R1)
     F.R1R2 <- q * (fR1R2/2 + fR1R1) + p * (fR1R2/2 + fR2R2)
     F.R2R2 <- q * (fR1R2/2 + fR2R2)
     F_seed[1, ] <- c(F.R1R1, F.R1R2, F.R2R2)
   }

   return(F_seed)
}



#' Infection Prevalence by Class
#' 
#' Calculate the infection prevalence grouped by class of a matrix 
#' of \code{\link{pine12}} solutions, and include total prevalence in the final
#' column. Time should be the rows and classes as columns.
#' 
#' Used internally during \code{\link{pine12}}.
#' 
#' @param x The population projection over time with solutions as rows and
#' classes as columns.
#' @param total Total prevalence. Must be calculated externally.
#' @return A vector of the relative prevalences of each of the classes, with
#' the total prevalence in the final column.
#' @author Stu Field
#' @seealso \code{\link{pine12}}
#' @examples
#' 
#' A <- matrix(1:96, ncol=12)
#' prevClass(A, total=1:8/10)
#' 
#' @export prevClass
prevClass <- function(x, total) {

   stages <- ncol(x)
   cols   <- stages / 2
   for ( i in 1:cols ) {
      if ( i==1 ) {
         prevMat <- x[, i + cols] / (x[, i] + x[, i + cols])
      } else {
         prevMat <- cbind(prevMat, x[, i + cols] / (x[, i] + x[, i + cols]) )
      }
   }
   dimnames(prevMat) <- list(paste0("Gen", 1:nrow(x)), paste0("Class", 1:cols))
   cbind(prevMat, total)
}


#' Group by Class (pine12)
#' 
#' A simple post-process function to group the susceptible & infected classes
#' in the 12-class pine12 model, resulting in 6 classes with solutions for each
#' time step of the projection. Summed by class gives: Seeds, SD1, SD2, SP, YA,
#' MA.
#' 
#' Function for calculating \code{\link{pine12}} stages by grouping infected 
#' and susceptible classes
#' 
#' @param x A storage matrix of the 12 class matrix projection containing
#' classes as cols and time steps as rows. Does not explicitly have to be 12
#' classes (cols) but for proper grouping, it must be only a SI model, i.e. two
#' disease classes only.
#' @return A matrix of the solutions over time grouped by stage (1-6).
#' @author Stu Field
#' @seealso \code{\link{pine12}}
#' @examples
#' 
#' # 8 generations/12 stages
#' A <- matrix(sample(1:96, 96, replace=TRUE), ncol=12) 
#' groupClass(A)
#' 
#' @export groupClass
groupClass <- function(x) {

   stages <- ncol(x)
   cols   <- stages / 2
   for ( i in 1:cols ) {
      if ( i==1 ) {
         groupMat <- x[, i] + x[, i + cols]
      } else {
         groupMat <- cbind(groupMat, x[, i] + x[, i + cols])
      }
   }
   dimnames(groupMat) <- list(paste0("Gen", 1:nrow(x)), paste0("Class", 1:cols))
   return(groupMat)
}



#' Vital Rates Calculation
#' 
#' Calculates survivorship and transition probabilities for a base model
#' classes given a vector of residence times and mortalities for each class.
#' 
#' @param R A vector of residence times for each class
#' @param M A vector of mortalities for each class. This is converted to
#' survivorship via 1 - m_i
#' @return A matrix with classes as the columns. Row 1 is the survivorship, row
#' 2 are transitions.
#' @note Used internally in \code{\link{pine12}}and \code{\link{pine36}}
#' @author Stu Field
#' @seealso \code{\link{pine12}}, \code{\link{pine36}}
#' @examples
#' 
#' M <- c(1, 0.152, 0.105, 0.02, 0.015, 0.005)
#' R <- c(1, 4, 16, 20, 50, 0)
#' vitals(R, M)
#' 
#' @export vitals
vitals <- function(R, M) {

   stopifnot( length(R) == length(M) )
   L <- length(R)
   gamma <- 1 / R
   sigma <- 1 - M
   surv  <- (1 - gamma) * sigma
   trans <- gamma * sigma
   surv[L]  <- sigma[L]  # residence time is infinite for adults
   trans[L] <- 0        # no transition rate for adults
   ret <- rbind(surv, trans)
   colnames(ret) <- paste0("Stage", 1:L)
   ret

}



#' Allelic Frequencies
#' 
#' Calculate allele frequencies (p, q) from population count
#' data. Number of individuals must be in the form c(AA, Aa, aa).
#' 
#' @param x Vector containing the counts for each of the three
#' genotypes: AA, Aa, aa
#' @return A vector of length 2 in the form c(p, q) or c(freq.A, freq.a).
#' @note further notes
#' @author Stu Field
#' @references Hardy-Weinberg and many others...
#' @examples
#' 
#' pop <- c(AA=123, Aa=45, aa=88)
#' freqz(pop)
#' 
#' @export freqz
freqz <- function(x) {
   A <- unname((2*(x[1]) + x[2]) / (sum(x)*2))
   a <- unname((2*(x[3]) + x[2]) / (sum(x)*2))
   c(A=A, a=a)
}



