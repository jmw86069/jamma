
#' Update function default parameters
#'
#' Update function default parameters
#'
#' This function is a minor extension to `update_list_elements()`
#' intended to help update function parameters which are defined
#' as a nested list. See examples.
#'
#' The main utility is for a function that defines a full set of
#' required argument values, where the user calling the function
#' may want to modify onyl a subset of those default values.
#'
#' @family jam utilities
#'
#' @param function_name `function` or `character` string referring to
#'    the name of a function, passed to `formals()`.
#' @param param_name `character` string of the argument/parameter name
#'    in `formals()` of the `function_name`.
#' @param new_values `list` or named vector used to update or augment
#'    existing default argument values defined for `param_name` of
#'    `formals()` of `function_name`.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' # function uses y as passed to the function
#' test_function_1 <- function(y=list(a=1, b=4)) {
#'    print("y:");
#'    print(y);
#' }
#' test_function_1(y=list(b=5, d=1:5))
#'
#' # function starts with y formals, adds or updates new values
#' test_function_2 <- function(y=list(a=1, b=4)) {
#'    y <- update_function_params(test_function,
#'       param_name="y",
#'       new_values=y);
#'    print("y:");
#'    print(y);
#' }
#' test_function_2(y=list(b=5:6, d=1:5))
#'
#' @export
update_function_params <- function
(function_name=NULL,
 param_name=NULL,
 new_values=NULL,
 verbose=FALSE,
 ...)
{
   ## Purpose is to facilitate updating default parameters which are present in a list,
   ## but where defaults are defined in a function name.  So if someone wants to change
   ## one of the default values, but keep the rest of the list of defaults, this function
   ## does it.
   ##
   ## The default values are taken from the function formals, using eval(formals(functionName)).
   default_params <- eval(formals(function_name)[[param_name]]);
   if (verbose) {
      jamba::printDebug("update_function_params(): ",
         "str(default_params)");
      print(str(default_params));
      jamba::printDebug("update_function_params(): ",
         "str(new_values)");
      print(str(new_values));
   }
   if (length(new_values) > 0) {
      default_params_list <- list(default_params);
      names(default_params_list) <- param_name;
      new_values_list <- list(new_values);
      names(new_values_list) <- param_name;
      default_params <- utils::modifyList(default_params_list,
         new_values_list)[[param_name]];
   }
   # default_params <- update_list_elements(default_params,
   #    new_values,
   #    verbose=verbose);
   return(default_params);
}

#' Update a subset of list elements (deprecated)
#'
#' Update a subset of list elements (deprecated)
#'
#' Deprecated.
#'
#' This function is intended to help update a nested `source_list`,
#' a subset of whose values should be replaced with entries
#' in `update_list`, leaving any original entries in `source_list`
#' which were not defined in `update_list`.
#'
#' This function may be useful when manipulating lattice or ggplot2
#' graphical parameters, which are often stored in a nested
#' list structure.
#'
#' Note: This function is deprecated, because there exists a function
#' `utils::modifyList()` written by Lattice co-author
#' Dr. Deepayan Sarkar, that serves the same purpose as this function.
#' No doubt the comment above is accurate, that this function is
#' widely used in the Lattice R package to update graphical
#' parameters.
#'
#' @family jam utilities
#'
#' @param source_list `list` with input values
#' @param update_list `list` with values to populate into `source_list`.
#' @param list_layer_num `integer` used for internal use by this function,
#'    which defines the list depth when traversing a nested list.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' # function starts with y formals, adds or updates new values
#' y <- list(a=1, b=4, e=NULL)
#' update_list_elements(y, list(b=5:6, d=1:5));
#'
#' # also see the preferred function utils::modifyList below
#' utils::modifyList(y, list(b=5:6, d=1:5));
#'
#' @export
update_list_elements <- function
(source_list,
 update_list,
 list_layer_num=1,
 verbose=FALSE,
 ...)
{
   ## Purpose is to update elements in a list, allowing for multi-layered lists
   ## of lists. In case of a list-of-list, it will call this function again with
   ## each successive layer of listedness.
   ##
   ## Handy for updating lattice graphics settings, which are impossibly nested
   ## tangled ball of textual yarn.  An example:
   ## tp2 <- updateListElements(trellis.par.get(), list(fontsize=list(points=6)));
   ## trellis.par.set(tp2);
   ##
   ## Or in one line:
   ## trellis.par.set(updateListElements(trellis.par.get(), list(fontsize=list(points=6))));
   if (length(update_list) == 0) {
      return(source_list);
   }
   if (class(update_list) %in% c("list") && class(update_list[[1]]) %in% c("list")) {
      for (update_list_name in names(update_list)) {
         if (update_list_name %in% names(source_list)) {
            ## If the name already exists, we must update items within the list
            source_list[[update_list_name]] <- update_list_elements(
               source_list=source_list[[update_list_name]],
               update_list=update_list[[update_list_name]],
               list_layer_num=list_layer_num+1);
         } else {
            ## If the name does not already exist, we can simply add it.
            source_list[[update_list_name]] <- update_list[[update_list_name]];
         }
      }
   } else {
      if (!is.null(names(update_list))) {
         source_list[names(update_list)] <- update_list;
      } else {
         source_list <- update_list;
      }
   }
   return(source_list);
}
