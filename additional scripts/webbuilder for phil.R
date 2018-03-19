# WebBuilder example for Phil
# not using the cheddar package


# source the webbuilder functions
source("Useful WebBuilder functions.R")
# function to make adjacency matrix from pairs
# used at end of script
source("pairs_to_adj function.R")


# registry of all NZ feeding links
registry <- read.csv("complete webbuilder registry.csv")

# new taxa that you want to build a web for
# e.g. a taxa list from one of your sites
# important that this is in the form of a data frame with the column titled "node"
# if not called "node" then WebBuilder throws an error
nodes <- read_csv("nodes for example.csv")

# taxonomoy information
# this is a large data frame with all of the taxonomic information for (most) freshwater taxa in NZ
# important to note the columns "res.method" and "con.method"
# This is where you can tweak the resolution
# e.g. res.method = genus --> res.method = family
taxonomy <- read_csv("taxonomy.csv")

# paste "minimum." to res/con.method
# this is another WebBuilder specific requirement
# if no columns in the "nodes" object are titled "minimum.res.method" WebBuilder throws an error
# could also obviously change this in the taxonomy csv
colnames(taxonomy)[5:6] <- paste("minimum.",
                          colnames(taxonomy[c(5:6)]),
                          sep = "")

# join the taxonomy information to your node list
nodes <- merge(nodes,
               taxonomy,
               by.x = "node", 
               by.y = "name")

# infer feeding links using WebBuilder()
# make sure that all of the methods in "minimum.res/con.method" columns are in the "method = c(...)" argument, otherwise throws an error
# e.g. in Gray example code they have "class", but since I don't use that I don't need it in here 
links <- WebBuilder(nodes,
                    registry,
                    method = c("exact",
                               "genus",
                               "family",
                               "order"))

# make adjacency matrix (binary predation matrix)
a.mat <- pairs_to_adj(nodes, links)

