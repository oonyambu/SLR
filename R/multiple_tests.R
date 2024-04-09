
#' Computes the tests between multiple combinations
#'
#' @param formula an object of class formula of the form y~factors or
#' y~factors|groups where y refers
#' to the response variable(s), factors determine the different test combination
#' and groups define the groups on which to carry the tests.
#' @param data a dataframe containing the variables defined in the formula
#' @param FUN the function f(x,y,...) used to carry out the test. eg t.test,
#' cor.test. default is t.test
#' @param ... additional arguments to the function FUN above
#' @param response the resulting name that defines the different responses in
#' the formula
#' @param var_name logical indicating whether the factor name should be
#' included when generating the combinations
#' @param select the quantity needed from the test carried out eg. p.value etc
#' @param wide logical indicating whether to return a long table or a wide table
#' @return a tidy dataframe containing all the quantities from the mutliple
#' tests
#'
#' @export
#'
#' @examples
#' multiple_tests(Sepal.Length~Species, iris)
#' # compare the first row to
#' t.test(Sepal.Length~Species, iris, subset = Species %in% c('setosa', 'versicolor'))
#' multiple_tests(cbind(Sepal.Length, Sepal.Width)~Species, iris)
#' multiple_tests(Sepal.Length + Sepal.Width~Species, iris)
#' multiple_tests(.~Species, iris)
#' multiple_tests(hp~am+vs, mtcars)
#' multiple_tests(qsec + hp~am+vs, mtcars, var_name = TRUE)
#' multiple_tests(hp~am|vs, mtcars, var_name = TRUE)
#' multiple_tests(.~Species, iris,cor.test)

multiple_tests <- function (formula, data, FUN = "t.test", ..., response = "response", 
                                           select = NULL, var_name = FALSE, wide = FALSE) {
  UseMethod('multiple_tests')
}

#' @export
multiple_tests.formula <- function (formula, data, FUN = "t.test", ..., response = "response", 
          select = NULL, var_name = FALSE,wide = FALSE) {
  FUN <- match.fun(FUN)
  fn <- function(x, g) {
    f <- function(cmbs) {
      nms <- chartr(".", "_", paste0(names(cmbs), collapse = ":"))
      cbind(grp = nms, my_tidy(FUN(cmbs[[1]], cmbs[[2]], 
                                   ...)))
    }
    do.call(rbind, combn(split(x, g), 2, f, simplify = FALSE))
  }
  fm <- function(x) {
    between <- if (var_name) 
      Map(paste, args$between, x[args$between], sep = "")
    else x[args$between]
    res <- lapply(x[args$response], fn, between)
    array2DF(structure(res, .Dim = length(res), dimnames = setNames(list(names(res)), 
                                                                    response)))
  }
  args <- from_formula(formula, data)
  result <- by(data, args$groups, fm)
  v <- lapply(model.frame(args$groups, data), unique)
  if (any(grepl("WITHIN", names(v))))
    names(v) <- "WITHIN"
  res <- array2DF(structure(result, .Dim = lengths(v), dimnames = v))
  rownames(res) <- NULL
  res$WITHIN <- NULL
  res$data.name <- NULL
  all_nms <- names(res)
  id_vars <- intersect(c(response, all.vars(args$groups)), all_nms)
  if (!is.null(select)) all_nms <- match.arg(c(id_vars,'grp',select), all_nms, TRUE)
  res <- res[,all_nms]
  if(wide) reshape(res, v.names = setdiff(all_nms,c(id_vars, 'grp')), 
                   timevar = 'grp', idvar = id_vars, sep = '_',
                   direction = 'wide')
  else res
}

# multiple_tests <- function(formula, data, FUN = 't.test', ...,
#                            response = 'response',
#                            select = NULL, var_name = FALSE){
#   FUN <- match.fun(FUN)
#   fn <- function(x, g){
#     f <- function(cmbs){
#       nms <- chartr('.', '_', paste0(names(cmbs), collapse = ":"))
#       cbind(grp = nms, my_tidy(FUN(cmbs[[1]], cmbs[[2]], ...)))
#     }
#     do.call(rbind, combn(split(x, g), 2, f, simplify = FALSE))
#     
#   }
#   
#   fm <- function(x){
#     between <-  if(var_name)  Map(paste, args$between, x[args$between], sep='')
#     else x[args$between]
#     res <- lapply(x[args$response], fn, between)
#     
#     array2DF(structure(res, .Dim = length(res),
#                        dimnames = setNames(list(names(res)), response)))
#   }
#   args <- from_formula(formula, data)
#   result <- by(data, args$groups, fm)
#   
#   v <- lapply(model.frame(args$groups, data), unique)
#   if(grepl('WITHIN', names(v))) names(v) <- "WITHIN"
#   res <- array2DF(structure(result, .Dim = lengths(v), dimnames = v))
#   rownames(res) <- NULL
#   res$WITHIN <- NULL
#   res$data.name <- NULL
#   subset(res,select = if(is.null(select))names(res)else match.arg(select, names(res),TRUE))
# }

from_formula <- function(form, data){
  if(!is.data.frame(data))
    stop('data must be a dataframe')
  if(length(form)!=3)
    stop("you must specify both the RHS and LHS of the formula", call. =FALSE)

  var_nms <- all.names(form, FALSE)
  if(anyDuplicated(var_nms))
    stop("cannot use any variable more than once", call. = FALSE)
  if('.' %in% all.vars(form[[3]]))
    stop("wrong formula. cannot have `.` on the RHS", call. = FALSE)

  LHS <- all.vars(form[[2]])
  if("." %in% LHS)  {
    if(missing(data))
      stop('wrong use of `.` without data')
    LHS <- c(setdiff(names(data), var_nms), LHS[LHS!='.'])
  }


  if(length(form[[3]])==3 && form[[3]][[1]] == '|'){
    RHS_BETWEEN <- all.vars(form[[3]][[2]])
    RHS_WITHIN <- reformulate(all.vars(form[[3]][[3]]))
  }
  else {
    RHS_BETWEEN <- all.vars(form[[3]])
    RHS_WITHIN <- ~cbind(WITHIN = rep(1, nrow(data)))
  }
  list(response = LHS, between = RHS_BETWEEN,
       groups = RHS_WITHIN)
}


searchsorted <- function(x, vec){
  idx <- findInterval(x, vec, all.inside = TRUE)
  vals <- vec[idx]
  idx2 <- abs(vals - x)>2
  vals2 <- vec[idx[idx2] + 1]
  is.na(vals2) <- vals2 - x[idx2]>2
  replace(vals, idx2, vals2)
}

my_tidy <- function(x){
  if(is.atomic(x)) return(data.frame(values = t(x)))
  if(!is.list(x)) return(x)
  fn <- function(i,j){
    if(length(i) == 1)
      data.frame(t(setNames(i, j)))
    else{
      nms <- if(j == 'conf.int') c('conf.low', 'conf.high')
      else paste0(j, seq_along(i))
      data.frame(t(setNames(i, nms)))
    }
  }
  x <- Filter(Negate(is.null), x)
  do.call(cbind,  Map(fn, unname(x), names(x)))
}

