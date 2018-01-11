# "supplemental fish calculations.R"
# estimating N~M relationship using taxa average data


# invertebrate biomass and abundance
invert <- readRDS("estimated invert bodymass.RDS")

lm_obj <- lm(log10(no.m2)~site:log10(dw), data = invert)
mean(lm_obj$coefficients[-1])
# mean slope = -0.3015064
range(lm_obj$coefficients[-1])
# range = [-0.5011 -0.1226]