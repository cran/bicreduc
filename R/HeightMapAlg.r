HeightMapAlg <- function(R,B)
{
    storage.mode(R) <- "double"
    storage.mode(B) <- "integer"
    .Call("HeightMapAlgorithm", R, B, PACKAGE="bicreduc")
}
