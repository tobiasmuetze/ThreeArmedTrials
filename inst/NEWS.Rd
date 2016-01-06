\name{NEWS}
\title{News for Package \pkg{ThreeArmedTrials}}
\section{Changes in Version 0.3-1}{
\subsection{Functionality Changes}{
\itemize{
\item Changed parameter file for permutation simulation and edited the simulation
      script to incoporate the respective changes
\item Fixed bug in function \code{test_RET_normal} that degrees of
      freedom are not calculated for variance of zero
}
}
\subsection{Minor Changes}{
\itemize{
\item R-packages \code{parallel}, \code{doMC}, and \code{foreach} are
      suggest to run the simulations
}
}
}
\section{Changes in Version 0.3-0}{
\itemize{
\item Renamed permutation test functions.
\item Added feature in \code{simulate_test_RET_permutation} to include
      competitor in simulation is new.
\item Added Wald-type test for normally distributed endpoints with
      heterogeneous variance.
}
}
\section{Changes in Version 0.2-1}{
\itemize{
\item Added simulation study for permutation test for binary
      and poisson endpoints to /inst/.
}
}
\section{Changes in Version 0.2-0}{
\itemize{
\item Included studentized permutation test \code{permutation.test} and Monte-Carlo
      simulation for the permutation test \code{simulation.permutation.test}.
\item Adjusted documentation.
\item Added default packages to import.
 }
}
\section{Changes in Version 0.1-0}{
\itemize{
\item Created.
 }
}



