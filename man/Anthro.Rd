\name{Anthro}
\docType{data}
\alias{Anthro}
\title{Children anthropometric data}
\description{
The dataset comprises 985 observations and 8 variables on children aged 0-60 months from the Central-West region of Brazil. The data are a cleaned subset of the 2006 Brazilian National Demographic Health Survey of Women and Children (NDHS), a complex survey in which sampling units were selected in two stages: primary sampling units (PSUs) were census tracts, and secondary sampling units (SSUs) were households. Tract selection within each stratum was designed to ensure a minimum number of blood samples, given the prevalence of vitamin A deficiency in children. Outlier weights (> 30 kg) and rows with missing values have been removed.
}
\usage{data("Anthro")}
\format{
  The data frame has the following components:
  \itemize{
    \item \bold{wgt} : weight
    \item \bold{hgt} : height
    \item \bold{wgt_ind} : weight-for-age index
    \item \bold{hgt_ind} : height-for-age index
    \item \bold{ruc} : rural-urban classification (urban = 1 and rural = 2) 
    \item \bold{sex} : boy = 1 and girl = 0
    \item \bold{age} : age in months
    \item \bold{dweight} : sampling weights
  }
}
\source{
The 2006 Brazilian National Demographic Health Survey of Women and Children (NDHS) was published by the Brazilian Institute of Geography and Statistics.
}

\keyword{datasets}

