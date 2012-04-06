.First <-
function () 
{
    options(repos = c(CRAN = "http://cran.r-project.org/"), browserNLdisabled = TRUE, 
        deparse.max.lines = 2)
}
.ls.objects <-
function (pos = 1, pattern, order.by = "Size", decreasing = TRUE, 
    head = TRUE, n = 10) 
{
    napply <- function(names, fn) sapply(names, function(x) fn(get(x, 
        pos = pos)))
    names <- ls(pos = pos, pattern = pattern)
    obj.class <- napply(names, function(x) as.character(class(x))[1])
    obj.mode <- napply(names, mode)
    obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
    obj.size <- napply(names, object.size)/10^6
    obj.dim <- t(napply(names, function(x) as.numeric(dim(x))[1:2]))
    vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
    obj.dim[vec, 1] <- napply(names, length)[vec]
    out <- data.frame(obj.type, obj.size, obj.dim)
    names(out) <- c("Type", "Size", "Rows", "Columns")
    out <- out[order(out[[order.by]], decreasing = decreasing), 
        ]
    if (head) 
        out <- head(out, n)
    out
}
