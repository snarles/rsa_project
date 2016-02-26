####
##  Generate data with similar structure as gambling data
##  Currently, ROIs are uncorrelated with each other, and each roi has the
##  same number of components
####

library(pracma)

####
#  A function with a lot of options!
#  Parameters:
#   npca : number of dimensions per roi
#   mixture_param : set to 0 to null model, 1 for alternative model
#  Specify:
#   'nrepeats'
#     AND
#   'params' OR 'ncopes' and 'q' (generates params matrix)
#     AND
#   'coeffs_true' OR  'q', nrois', 'npca',
#    'nsubjects' (generates random params) and 'mixture_param'
#     AND
#   'Sigmas_true' OR  'nrois', 'npca', 'nsubjects' and 'Wishart_mult'
#     (uses different Sigmas for each subject)
generate_copes_data <- function(nrepeats = 3, 
                                nsubjects = NULL, ncopes = NULL, nrois = NULL,
                                npca = NULL, mixture_param = NULL,
                                Wishart_mult = NULL,
                                coeffs_true = NULL, Sigmas_true = NULL,
                                ) {
  
}