################################################################################
### Compute the intersection of a circular domain with a polygonal domain of
### various classes (currently: "owin" or "SpatialPolygons")
###
### Copyright (C) 2009-2015 Sebastian Meyer
###
### This file is part of the R package "surveillance",
### free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at https://www.R-project.org/Licenses/.
################################################################################


intersectPolyCircle.owin <- function (object, center, radius, npoly = 32, ...)
{
    circle <- disc(radius = radius, centre = center, npoly = npoly)
    res <- intersect.owin(circle, object)  # order does not affect runtime
    ## ensure "polygonal" type (because of rescue.rectangle in intersect.owin)
    as.polygonal(res)
}

intersectPolyCircle.SpatialPolygons <- function (object, center, radius,
                                                 npoly = 32, ...)
{
    circle <- discpoly(center, radius, npoly = npoly, class = "Polygon")
    circleSpP <- SpatialPolygons(list(Polygons(list(circle), "0")))
    ## ensure that circleSpP has exactly the same proj4string as 'object'
    circleSpP@proj4string <- object@proj4string
    rgeos::gIntersection(circleSpP, object)
}
