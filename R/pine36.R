##################################
### 36 CLASS GENETICS MODEL ######
### with Infection ###############
### with Genetics 36x36 ##########
### June 2nd 2011 ################
##################################

#' The pine36 Genetic Infection Model
#' 
#' The full 36 class model with susceptible and infected in the pine12 model. 
#' There are three forms of nonlinearity: 2 forms of density-dependent seed 
#' germination (rALs, r.cache) & density-dependent fecundity (r.cones). In 
#' the model, infection has a fitness cost to \emph{both} survivorship 
#' and fecundity.
#'
#' Genetics of resistance is included the model via simple dominance form 
#' of inheritance where the resistance allele R1 [p] is initially rare and 
#' fully dominant to the susceptible allele (R2; q). Fitness is modeled as 
#' partial dominance where the \code{h} argument is the coefficient 
#' specifying the degree of dominance (see arguments). A value of 0 (default) 
#' represents the scenario of complete
#' dominance where heterozygotes and homozygotes for the resistance allele
#' have identical fitness. When \code{h=1}, heterozygote fitness is identical
#' to the susceptible homozygotes, effectively negating resistance. When \code{h=0.5}
#' heterozygote fitness is intermediate between the two extremes.
#' 
#' @param Gen Integer. Number of generations to project the population
#' @param x1 A matrix indicating the initial population, with rows as classes
#' and cols as genotypes
#' @param M Numeric. A vector (of length=5) describing the mortality for each class (1
#' -> 5). Subsequently converted to survivorship by the \code{vitals}
#' functions
#' @param m6 The mortality of the adult class (6)
#' @param R Numeric. A vector (of length=6) describing the mean Residence time for each
#' class (1 -> 6)
#' @param Beta Numeric scalar. The transmission probability
#' @param dbh.v A vector (of length=6) describing the mean dbh for each class
#' (1 -> 6)
#' @param s1 Cost of infection to survivorship of SD1 individuals
#' @param s2 Cost of infection to survivorship of SD2 individuals
#' @param delta Effect of infection on survivorship
#' @param alpha1 Coefficient corresponding to the LAI of the sapling (SA)
#' stage
#' @param alpha2 Coefficient for the conversion of dbh -> LAI
#' @param alpha3 Coefficient in the exponent for the conversion of dbh -> LAI
#' @param LAIb Background LAI due to external species in the simulation
#' @param Cmax Maximum cone production per reproductive tree
#' @param S.cone Seeds per cone
#' @param P.find Proportion of seeds found by dispersal vectors (birds)
#' @param P.cons Proportion of seeds immediately consumed by birds
#' @param SpC Seeds per cache
#' @param cf Cost of infection to fecundity
#' @param nBirds Number of birds per modeling plot (area)
#' @param r.site Suitability of the site to germination. Defaults to no effect
#' (=1)
#' @param h Numeric. Scalar in [0,1] indicating the degree of dominance for the
#' resistance allele
#' @param rho Relative difference in fecundity between young adults (jv) &
#' mature adults (ad)
#' @param qpollen Logical or Numeric. If \code{FALSE} (default), a pollen cloud,
#' i.e. a \emph{constant} q allele, is turned not assumed and offspring are produced 
#' according to \code{mating36}. If Numeric, the frequency of the susceptible allele 
#' (q/R2) which is fixed. Must be in [0,1]
#' @param NoSeln Logical. Should there be no selection applied to both
#' survivorship & fecundity? Defaults to \code{FALSE} and therefore WITH
#' selection. Used for testing Hardy-Weinberg genetics
#' @param plot Logical. Should output plots be created?
#' @param csv Logical. Should all solutions and time steps be converted to a
#' 2-dim matrix and exported to a *.csv file (in the working directory)?
#' @param filename Character string. If \code{csv=TRUE}, the name of the csv
#' file to be produced
#' @param silent Logical. Should the output be silenced?
#' @return A list containing:
#' \item{theta }{List of model parameters}
#' \item{InitialPop }{A matrix corresponding to \code{x1}, the initial population}
#' \item{FinalPopn }{A matrix of the final solution of the projection}
#' \item{Gtypes }{A matrix of the solutions of all time steps, summed by Genotype}
#' \item{FinalSum }{The sum of the population at the final time step}
#' \item{PopTotals }{A vector of the total population by over each time step 
#' in the projection}
#' \item{R2 }{A vector of the frequency of the R2/q allele (susceptible) 
#' at each time step in the projection}
#' \item{GerminationRate }{A vector of the germination rate (density-dependent) at 
#' each time step in the projection}
#' \item{r.cones }{The value of \code{r.cones} at the final time step in
#' the projection}
#' \item{FinalSeeds }{A matrix of the number of seeds in each genotype class
#' added to the population in the final step of the projection. Non-zero
#' entries only in the first row}
#' \item{LAIy }{A vector of LAI of the intermediate population, following
#' survivorship & transition but prior to fecundity & infection}
#' \item{Prevalence }{A vector of the total infection prevalence 
#' (proportion infected) at each time step in the projection}
#' \item{LambdaGrowth }{A vector of the proportional growth in successional time
#' steps during the projection}
#' \item{LambdaGrowthGtype }{A matrix of the proportional growth in
#' successional time steps during the projection, separated by genotype}
#' \item{FullSolution }{An array of the population projection with 
#' rows = classes, cols = genotypes, planes = time steps}
#' @note Lots of work during my post-doc!
#' @author Stu Field
#' @seealso \code{\link{pine12}}
#' @references %% ~put references to the literature/web site here ~
#' @examples
#' 
#' R1R1 <- c(10000,100,87,82,123,244,   0,0,0,0,0,0)
#' #R1R1 <- c(0,0,0,0,0,0,   0,0,0,0,0,0)  # fresh colonization
#'
#' R1R2 <- c(10000,100,87,82,123,244,   0,0,0,0,0,0) * 2
#' #R1R2 <- c(0,0,0,0,0,0,   0,0,0,0,0,0)
#'
#' R2R2 <- c(10000,100,87,82,123,244,   0,0,0,0,0,0)
#' #R2R2 <- c(0,0,0,0,0,0,   0,0,0,0,0,0)
#'
#' pine36(Gen=25, filename="Run1", plot=TRUE, silent=FALSE)
#' 
#' @importFrom stuRpkg diagR zeros BlockMat
#' @export pine36
pine36 <- function(Gen = 25,
                   x1 = cbind(c(10000,100,87,82,123,244, 0,0,0,0,0,0),
                              c(10000,100,87,82,123,244, 0,0,0,0,0,0) * 2,
                              c(10000,100,87,82,123,244, 0,0,0,0,0,0)),
                   M = c(1, 0.152, 0.105, 0.02, 0.015),
                   m6 = 0.005,
                   R = c(1, 4, 16, 20, 50, 0), 
                   Beta = 0.044,
                   dbh.v = c(0, 0, 0, 2.05, 12.5, 37),
                   #s1 = 0.01, s2= 0.13,
                   s1 = 0.99,
                   s2 = 0.87,
                   delta = 0.15,
                   alpha1 = 0.456,
                   alpha2 = 0.0736,
                   alpha3 = 2.070,
                   LAIb = 0,
                   Cmax = 7.5,
                   S.cone = 46,
                   P.find = 0.8,
                   P.cons = 0.3,
                   SpC = 3.7,
                   cf = 0.875,
                   nBirds = 3,
                   r.site = 1,
                   h = 0,
                   rho = 0.1,
                   qpollen = FALSE,
                   NoSeln = FALSE,
                   plot = FALSE,
                   csv = FALSE,
                   filename=NULL,
                   silent=TRUE) {
   #########################################
   # Storage array and  vectors for results
   #########################################
   n_stages <- length(x1)
   classes  <- c(paste0("C",1:6), paste0("C",1:6,"i"))
   I.Gtype  <- apply(x1[-1,], 2, sum)

   # overall storage of population projection
   Popn <- array(0, dim=c(n_stages/3, 3, Gen),
                 dimnames=list(classes, names(I.Gtype), NULL))

   # first plane is the initial population
   Popn[,,1] <- x1              

   # storage of intermediate popn following survival 
   # used only for debugging; not returned
   Popn.S <- array(0, dim=dim(Popn), dimnames=dimnames(Popn))

   # Setup storage vectors for variables of interest ###
   freq.R         <- freqz(I.Gtype)[1]
   freq.r         <- freqz(I.Gtype)[2]
   Gtype.Sum      <- matrix(NA, ncol=3, nrow=Gen, dimnames=list(NULL,names(I.Gtype)))
   Gtype.Sum[1, ] <- I.Gtype                           # sum pop by genotype
   Pop.Total      <- sum(x1[-1, ])                     # Initial PopSize (remove seeds)
   Prev.v         <- sum(x1[8:12, ]) / (sum(x1[-1, ])) # Initial prevalence
   SR.v           <- numeric(Gen)                      # storage vector Seedling Recruitment
   LAI.v          <- numeric(Gen)                      # storage vector LAI
   r.cache.v      <- numeric(Gen)                      # storage vector r.cache


   ########################################
   # Calculate Vital Rates (vitals function)
   ########################################
   M <- c(M, m6)
   S.par <- vitals(R, M)["surv", ]
   T.par <- vitals(R, M)["trans", ]

   #######################
   ## Set linear Parameters
   #######################
   
   # Leaf Area Calculation
   LA.v <- alpha2 * (dbh.v^alpha3)
   LA.v[3]= alpha1               ### Don't forget the sedondary seedlings

   # Fitness pars
   # w[1] = R1R1
   # w[2] = R1R2
   # w[3] = R2R2
   # Viability Cost - Survivorship & Transition
   #s.vec <- c(1, s1, s2, exp(-delta*(dbh.v[4:6])) )
   s.vec <- c(0, s1, s2, exp(-delta*(dbh.v[4:6])) ) # Doesn't matter which; multiplied by zero because seeds don't survive the matrix multiplication.

   I6   <- diag(6)
   w.hs <- diag(1 - s.vec * h)
   w.s  <- diag(1 - s.vec)

   # Beta parameters
   B.vec <- c(0, rep(Beta,5))    # Beta is constant for all classes

   # Set up projection matrix (Linear Map)
   # Survivorship Matrix
   S <- diag(S.par) + diagR(T.par[1:5], -1)

   # Fitness Matrix
   W.AA <- BlockMat(list(I6,zeros(6), zeros(6),I6), b=2)
   W.Aa <- BlockMat(list(I6,zeros(6), zeros(6),w.hs), b=2)
   W.aa <- BlockMat(list(I6,zeros(6), zeros(6),w.s), b=2)

   # Zeros Matrix for filling
   Z <- zeros(n_stages / 3)

   blocks <- list( W.AA, Z,    Z,
                   Z,    W.Aa, Z,
                   Z,    Z,    W.aa)
               
   W.s <- BlockMat(blocks, b=3)

   if ( NoSeln )    # No selection
      W.s <- diag(36)

   # Beta Matrix
   B <- diag(B.vec)

   # Combine sub-matrices for survival
   SM12 <- BlockMat(list(S, zeros(6), zeros(6), S), b=2)
   SM36 <- BlockMat(list(SM12,Z,Z, Z,SM12,Z, Z,Z,SM12), b=3)

   # Combine sub-matrices for infection
   BM12 <- BlockMat(list(I6-B,zeros(6), B,I6), b=2)
   BM36 <- BlockMat(list(BM12,Z,Z, Z,BM12,Z, Z,Z,BM12), b=3)

   # Loop over generations
   for ( n in 2:Gen ) {

      # Reset x.vec for new iteration
      x.vec <- Popn[,,n-1]
      
      # LAI.x Calculation
      # Projected LAI of last year's popn (x)
      # 1 ha = 10000 m^2
      LAI.x <- sum(LA.v * as.vector(x.vec)) / 10000 + LAIb

      # Determine New Seedlings
      SpB <- sum(x.vec[1,]) / nBirds
      r.cache <- 0.73 / (1 + exp((31000 - SpB)/3000)) + 0.27
      r.cache.v[n] <- r.cache               # Store r.cache
      r.ALs <- 1 / (1 + exp(2*(LAI.x - 3)))

      # Seedling transition (proportion seeds becoming seedlings)
      SR <- (((1 - P.find) * (1 - P.cons)) / SpC) * r.cache * r.ALs * r.site
      SR.v[n] <- SR                        # Store SR (seedline recruitment)
      New.Seedlings <- Popn[1,,n-1] * SR   # Determine number of seedlings

      # Add seedlings to population
      x.vec[2,] <- x.vec[2,] + New.Seedlings

      # Survive & transition
      # POST-multiplication of cost matrix W.s
      y.vec <- (SM36 %*% W.s) %*% as.vector(x.vec)
      
      Popn.S[,,n] <- y.vec     # Matrix stores survivorship vectors at each iteration

      # LAI.y Calculation
      LAI.y    <- sum(LA.v * y.vec) / 10000 + LAIb   # LAI of this year's popn
      LAI.v[n] <- LAI.y                              # Store LAI over time
      
      # Determine New Seeds
      r.cones <- (0.5/(1 + exp(5*(LAI.y - 2.25))) + 0.5)
      C.tree  <- r.cones * Cmax

      # Set up matrix of Fecundities
      # Fecundity Matrix
      F.mat <- matrix(0, nrow=n_stages/3, ncol=3)
      F.mat[c(5,11),] <- rho * (S.cone * C.tree)
      F.mat[c(6,12),] <- S.cone * C.tree

      # Fitness Matrix
      W <- matrix(0, nrow=n_stages/3, ncol=3)
      W[5:6,]   <- 1                                    # no cost to uninfected genotypes
      W[11:12,] <- rep(c(1, 1 - h * cf, 1-cf), each=2)  # fitness of infected gtypes
      if ( NoSeln )
         W[11:12,] <- 1                   # No selection
      
      # Hadamard product
      #W.f <- W * F.mat
      
      # convenient matrix Y of the surviving population; cols = genotypes
      Y <- matrix(y.vec, ncol=3)
      
      #######################
      ###     Mating      ###
      #######################
      # This uses Jesse's simplified code from FEScUE
      # assumes explicitely random mating
      # outputs matrix with next year's seeds
      ###########################################
      if ( !qpollen )
         Seeds <- mating36(Y, W, F.mat)
      if ( qpollen > 0 )
         Seeds <- mating36p(Y, W, F.mat, qpollen) # fixed pollen cloud
      
      # Add newborns to population
      y.vec <- as.vector(Y + Seeds)      # Add seeds to population, in the Fall (timing is important!)

      # Infect the population
      Popn[,,n] <- BM36 %*% y.vec          # Infect in Fall (last operation)


      ################
      #   OUTPUTS    #
      ################
      Gtypes <- apply(Popn[-1,,n], 2, sum)  # sum the cols (Gtypes); 
                                            # put [-1,x,n] to remove seeds
      Gtype.Sum[n,] <- Gtypes
      freq.R[n] <- freqz(Gtypes)[1]         # Calc p
      freq.r[n] <- freqz(Gtypes)[2]         # Calc q
      Pop.Total[n] <- sum(Popn[-1,,n])      # remove seeds
      Prev.v[n] <- ifelse(h==0,
                        sum(Popn[8:12,3,n]) / (sum(Popn[-1,,n])), 
                        sum(Popn[8:12,2:3,n]) / (sum(Popn[-1,,n])))
   }
   # END OF LOOP

   # Model parameters
   theta <- list(Gen = Gen,
                 Transmission = Beta,
                 LAIb = LAIb, 
                 r.site = r.site,
                 delta = delta, 
                 cf = cf,
                 Cmax = Cmax, 
                 SpC = SpC, 
                 S.cone = S.cone,
                 P.find = P.find, 
                 nBirds = nBirds, 
                 P.cons = P.cons,
                 alpha1 = alpha1, 
                 alpha2 = alpha2,
                 alpha3 = alpha3,
                 s1 = s1, 
                 s2 = s2,
                 h = h,
                 rho = rho,
                 qpollen = qpollen,
                 NoSeln = NoSeln,
                 plot = plot,
                 csv = csv,
                 silent = silent
            )

   theta$tree_pars <- data.frame(rbind(dbh.v=dbh.v, M=M, R=R))
   names(theta$tree_pars) <- paste0("Stage", seq_along(dbh.v))

   
   #############################
   # Post-process variables
   #############################
   # Post-process Lambda (Proportional Growth)
   Lambda  <- LambdaGrow(Pop.Total)              # Total population Lambda
   Lambda2 <- apply(Gtype.Sum, 2, LambdaGrow)    # Lambda by genotype


   #########################################
   # Convert Solutions & Export to CSV
   #########################################
   if ( csv ) {
      if ( is.null(filename) )
         stop("Please pass a file name for the CSV output file via the filename= argument")
      full2x2 <- do.call(rbind,
                         lapply(1:Gen, function(i) rbind(as.vector(Popn[,,i]))))
      rownames(full2x2) <- paste0("t", 1:Gen)
      colnames(full2x2) <- paste0("Class_", 1:n_stages)
      write.csv(full2x2, file=sprintf("%s_FullSoln_%s.csv", filename, Sys.Date()))
   }


   ############################
   # Plot variables of interest
   # one day pull this plotting routine out into separate function
   ############################
   if ( plot ) {
      par(mfrow=c(2, 3))
      plot(1:Gen, Pop.Total, type="l", 
           xlim=c(0,Gen), 
           xlab="Years", 
           main="Total Population", 
           ylab="Total Whitebark Population (no seeds)", 
           col="navy", lwd=2)
      grid()
      plot(1:Gen, freq.r, 
           xlab="Years", 
           main="R2 Allele", 
           ylab="q", 
           col="darkgreen")
      grid()
      plot(1:Gen, Prev.v, type="l", 
           col="darkred", 
           lwd=2, 
           xlab="Years", 
           main="Rust Prevalence", 
           ylab="Proportion infected individuals")
      grid()
      plot(1:(Gen-1), LAI.v[2:Gen], 
           xlab="Years", 
           main="Post-Survival Leaf Area Index", 
           ylab="LAI.y")
      grid()
      plot(1:(Gen-1), SR.v[2:Gen], 
           xlab="Years", 
           main="Germination Rate", 
           ylab="SR")
      grid()
      plot(1:(Gen-1), r.cache.v[2:Gen], 
           xlab="Years", 
           main="r.cache", 
           ylab="r.cache")
      grid()
   }

   ret <- list(theta = theta,                # model parameters
               InitialPop = x1,              # initial popn
               FinalPopn = Popn[,,Gen],      # Final Solution
               Gtypes = Gtype.Sum,           # Solutions summed by Genotype
               FinalSum = Pop.Total[Gen],    # Total Popn
               PopTotals = Pop.Total,        # Total Popn vector over time
               R2 = unname(freq.r),          # Frequency of 'q' allele final time (t)
               GerminationRate = SR.v,       # Vector of Germination rates
               r.cones = r.cones,            # r.cones in final time (t)
               FinalSeeds = Seeds,           # Seeds produced in final time (t)
               LAIy = LAI.v,                 # LAI.y in final time (t)
               Prevalence = Prev.v,          # Vector of disease prevalence
               LambdaGrowth = Lambda,        # Lambda = x(t) / x(t-1)
               LambdaGrowthGtype = Lambda2,  # Lambda by Genotype
               FullSolution = Popn)          # ALL Solutions by Generation

   if ( silent )
      invisible(ret)
   else
      ret

}

