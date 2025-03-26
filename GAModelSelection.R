# =============================================================================
# GAModelSelection Package Functions
# =============================================================================
# This script implements core functions for genetic algorithm–based model selection:
#   - preprocess_data(): Handles missing values, normalization, and optional shifting.
#   - generate_shift_ranges(): Generates a sequence of shift values.
#   - adjust_dep_var_for_family(): Adjusts the dependent variable to meet distribution requirements.
#   - build_model(): Constructs a candidate model based on a configuration.
#   - evaluate_fitness(): Computes fitness using DHARMa diagnostics.
#       * If all four diagnostic p-values exceed (0.05 - relaxation_term), are non-NA, and their product is nonzero,
#         then Fitness = 10000 + (product of p-values) + logLik(model); otherwise, Fitness = -Inf.
#   - run_model_selection(): Runs candidate generation, model building, diagnostics, and selection.
#       It compiles a diagnostic table that now includes information about all genes, discards candidates with bad fitness,
#       and returns a list containing:
#         * all_results: Full candidate evaluations.
#         * diagnostic_table: Data.frame summarizing candidates that pass the threshold.
#         * best_model: The model from the candidate with the highest fitness.
#         * best_model_sim: The DHARMa simulation object for the best model.
#         * top_models: List of top candidate evaluations (up to 10) that passed the threshold.
#         * passed_models: List of candidate evaluations that passed the fitness threshold.
#         * model_list: The full list of candidate evaluations (indexed by candidate ID).
#   - demo_model_selection(): Demonstrates the full workflow using the mtcars dataset.
#
# Note: For mixed-model examples with proper random effects, consider using datasets like lme4's sleepstudy.
# =============================================================================

# Load required packages
library(lme4)
library(glmmTMB)
library(DHARMa)
library(bestNormalize)
library(MuMIn)
library(boot)
library(caret)
library(car)
library(pbapply)  # For progress bar

#' Preprocess Data
#'
#' Preprocesses the dataset by handling missing values, applying normalization, and optionally adding a small shift.
#'
#' @param data A data.frame.
#' @param dep_var Character. Name of the dependent variable.
#' @param normalize_method Character. One of "no_transform", "arcsinh_x", "boxcox", "log_x", 
#'        "double_reverse_log", "sqrt_x", "yeojohnson", "Lambert".
#' @param missing_values_method Character. One of "none", "impute_mean", "discard".
#' @param shift_addition Logical. Whether to add a small shift.
#' @param step_for_shift Numeric. Step size for generating shift ranges.
#'
#' @return A preprocessed data.frame.
#' @export
preprocess_data <- function(data, dep_var, normalize_method = "no_transform",
                            missing_values_method = "none", shift_addition = FALSE,
                            step_for_shift = 0.01) {
  if (missing_values_method == "impute_mean") {
    for (col in names(data)) {
      if (is.numeric(data[[col]])) {
        data[[col]][is.na(data[[col]])] <- mean(data[[col]], na.rm = TRUE)
      }
    }
  } else if (missing_values_method == "discard") {
    data <- na.omit(data)
  }
  
  if (normalize_method != "no_transform") {
    norm_obj <- bestNormalize::bestNormalize(x = data[[dep_var]], allow_lambert = TRUE)
    data[[dep_var]] <- predict(norm_obj)
  }
  
  if (shift_addition) {
    data[[dep_var]] <- data[[dep_var]] + 1e-5
  }
  
  return(data)
}

#' Generate Shift Ranges
#'
#' Generates a numeric vector of shift values from 1e-5 to 1.
#'
#' @param step_for_shift Numeric. Step size (default 0.01).
#'
#' @return A numeric vector.
#' @export
generate_shift_ranges <- function(step_for_shift = 0.01) {
  seq(1e-5, 1, by = step_for_shift)
}

#' Adjust Dependent Variable for Selected Family
#'
#' Checks and adjusts the dependent variable to meet the requirements of the specified family.
#'
#' @param data A data.frame.
#' @param dep_var Character. Name of the dependent variable.
#' @param family Character. The selected family (e.g., "poisson", "Gamma", "beta_family").
#'
#' @return A data.frame with the adjusted dependent variable.
#' @export
adjust_dep_var_for_family <- function(data, dep_var, family) {
  y <- data[[dep_var]]
  family_lower <- tolower(family)
  
  if (family_lower == "poisson") {
    if (any(y < 0)) {
      shift_value <- abs(min(y)) + 1
      message("Shifting dependent variable by ", shift_value, " for Poisson family.")
      data[[dep_var]] <- y + shift_value
    }
    data[[dep_var]] <- round(data[[dep_var]])
  } else if (family_lower == "gamma") {
    if (any(y <= 0)) {
      shift_value <- abs(min(y)) + .Machine$double.eps
      message("Shifting dependent variable by ", shift_value, " for Gamma family.")
      data[[dep_var]] <- y + shift_value
    }
  } else if (family_lower %in% c("beta_family", "beta")) {
    if (any(y < 0) || any(y > 1)) {
      message("Scaling dependent variable for Beta family: transforming values to (0,1) range.")
      y_min <- min(y)
      y_max <- max(y)
      if (y_max - y_min == 0) stop("Cannot scale constant dependent variable for Beta family.")
      data[[dep_var]] <- (y - y_min) / (y_max - y_min)
    }
  }
  return(data)
}

#' Build Model from Candidate Genotype
#'
#' Constructs a statistical model based on the candidate's configuration, a fixed-effects formula,
#' and (optionally) random effects.
#'
#' @param candidate A list with candidate settings.
#' @param data A data.frame.
#' @param formula_base A formula for fixed effects.
#' @param random_effects Character vector for random effects.
#'
#' @return A fitted model object, or NULL on error.
#' @export
build_model <- function(candidate, data, formula_base, random_effects = NULL) {
  model_approach <- candidate$model_approach
  family_str <- candidate$family
  
  if (family_str %in% c("gaussian", "binomial", "poisson", "Gamma")) {
    family_obj <- get(family_str, mode = "function")()
  } else {
    family_obj <- family_str
  }
  
  fixed_formula <- formula_base
  if (!is.null(random_effects) && length(random_effects) > 0) {
    re_terms <- paste0("(1|", random_effects, ")")
    full_formula <- as.formula(paste(deparse(fixed_formula), paste(re_terms, collapse = " + "), sep = " + "))
  } else {
    full_formula <- fixed_formula
  }
  
  model <- tryCatch({
    if (model_approach == "lmer") {
      lme4::lmer(full_formula, data = data)
    } else if (model_approach == "glmer") {
      lme4::glmer(full_formula, data = data, family = family_obj)
    } else if (model_approach == "glmer.nb") {
      lme4::glmer.nb(full_formula, data = data)
    } else if (model_approach == "glmmTMB") {
      glmmTMB::glmmTMB(full_formula, data = data, family = family_obj)
    }
  }, error = function(e) {
    message("Error in build_model: ", e$message)
    NULL
  })
  return(model)
}

#' Evaluate Model Fitness
#'
#' Evaluates a fitted model using DHARMa diagnostics and computes a fitness score.
#' If all four diagnostic p-values (KS, dispersion, outlier, hetero) exceed (0.05 - relaxation_term),
#' and none are NA and their product is nonzero, then:
#'   Fitness = 10000 + (product of p-values) + logLik(model)
#' Otherwise, Fitness = -Inf.
#'
#' @param model A fitted model object.
#' @param relaxation_term Numeric. Adjustment for the p-value threshold.
#'
#' @return A list with elements:
#'   \item{fitness}{Numeric fitness value.}
#'   \item{p_values}{Named vector of diagnostic p-values.}
#'   \item{prod_p}{Product of the 4 p-values.}
#'   \item{sim}{DHARMa simulation object.}
#' @export
evaluate_fitness <- function(model, relaxation_term = 0) {
  sim_obj <- DHARMa::simulateResiduals(fittedModel = model, plot = FALSE)
  p_ks <- DHARMa::testUniformity(sim_obj, plot = FALSE)$p.value
  p_dispersion <- DHARMa::testDispersion(sim_obj, plot = FALSE)$p.value
  p_outlier <- DHARMa::testOutliers(sim_obj, plot = FALSE)$p.value
  p_hetero <- DHARMa::testGeneric(
    sim = sim_obj,
    summary = function(residuals) {
      groups <- cut(sim_obj$fittedPredictedResponse, breaks = 5)
      lev_result <- suppressWarnings(suppressMessages(car::leveneTest(residuals ~ groups)))
      lev_result[["Pr(>F)"]][1]
    },
    alternative = "two.sided",
    methodName = "Levene's Test for Heteroscedasticity",
    plot = FALSE
  )$p.value
  
  p_vals <- c(ks = p_ks, dispersion = p_dispersion, outlier = p_outlier, hetero = p_hetero)
  prod_p <- p_ks * p_dispersion * p_outlier * p_hetero
  logLik_value <- as.numeric(logLik(model))
  
  if (any(is.na(p_vals)) || prod_p == 0) {
    fitness <- -Inf
  } else if (all(p_vals > (0.05 - relaxation_term))) {
    fitness <- 10000 + prod_p + logLik_value
  } else {
    fitness <- -Inf
  }
  
  return(list(fitness = fitness,
              p_values = p_vals,
              prod_p = prod_p,
              sim = sim_obj))
}

#' Run Model Selection Using a Genetic Algorithm Approach
#'
#' Orchestrates the model selection process: preprocesses data, generates candidate configurations,
#' adjusts the dependent variable per candidate, builds candidate models, evaluates their fitness, and compiles a
#' diagnostic table. Candidates with fitness = -Inf, NA p-values, or a zero product are discarded.
#'
#' The diagnostic table includes the following columns:
#'   - Candidate: Candidate ID.
#'   - ModelApproach: Candidate's model approach.
#'   - Family: Candidate's family.
#'   - OutlierRemoval: Candidate's outlier removal flag.
#'   - RangesOutlierRemoval: Candidate's range for outlier removal.
#'   - OutlierFormula: Candidate's outlier formula.
#'   - ZI_Formula: Candidate's zero-inflation formula.
#'   - NormalizeApproach: Candidate's normalization approach.
#'   - ShiftAddition: Candidate's shift addition flag.
#'   - Link: Candidate's link function.
#'   - OutlierRemovalApproach: Candidate's outlier removal approach.
#'   - ZeroInflationApproach: Candidate's zero inflation approach.
#'   - Regularization: Candidate's regularization option.
#'   - FeatureSelection: Candidate's feature selection option.
#'   - GA_StoppingCriterion: Candidate's GA stopping criterion.
#'   - PopulationInitialization: Candidate's population initialization option.
#'   - CV_Folds: Candidate's cross-validation folds.
#'   - Fitness: Computed fitness value.
#'   - ProdP: Product of the 4 p-values.
#'   - KS, Dispersion, Outlier, Hetero: The four diagnostic p-values.
#'   - OutlierThreshold: The candidate's outlier threshold (from ranges_of_outlier_removal).
#'
#' Additionally, the function returns a list with:
#'   - all_results: Full candidate evaluations.
#'   - diagnostic_table: Data.frame of candidates that pass the threshold.
#'   - best_model: The model from the candidate with the highest fitness.
#'   - best_model_sim: DHARMa simulation object for the best model.
#'   - top_models: List of top candidate evaluations (up to 10).
#'   - passed_models: List of candidate evaluations that passed the fitness threshold.
#'   - model_list: The full list of candidate evaluations.
#'
#' @param data A data.frame.
#' @param dep_var Character. Dependent variable name.
#' @param predictors Character vector of predictor names.
#' @param random_effects Character vector for random effects.
#' @param formula_base A formula for fixed effects (if NULL, a default is constructed).
#' @param genetic_map A list defining available genes and alleles.
#' @param GA_control List of GA parameters (popSize, maxiter, run, parallel).
#' @param relaxation_term Numeric. Adjustment for p-value thresholds.
#' @param step_for_shift Numeric. Step size for generating shift ranges.
#' @param missing_values_method Character. One of "none", "impute_mean", "discard".
#' @param normalize_method Character. Normalization method.
#'
#' @return A list containing:
#'   - all_results: List of candidate evaluations.
#'   - diagnostic_table: Data.frame summarizing candidates that pass the threshold.
#'   - best_model: The model from the candidate with the highest fitness.
#'   - best_model_sim: The DHARMa simulation object for the best model.
#'   - top_models: List of top candidate evaluations (up to 10).
#'   - passed_models: List of candidate evaluations that passed the fitness threshold.
#'   - model_list: The full list of candidate evaluations.
#' @export
run_model_selection <- function(data, dep_var, predictors, random_effects = NULL,
                                formula_base = NULL, genetic_map,
                                GA_control = list(popSize = 50, maxiter = 100, run = 50, parallel = TRUE),
                                relaxation_term = 0, step_for_shift = 0.01,
                                missing_values_method = "none", normalize_method = "no_transform") {
  data <- preprocess_data(data, dep_var, normalize_method = normalize_method,
                          missing_values_method = missing_values_method,
                          shift_addition = genetic_map$shift_addition[1],
                          step_for_shift = step_for_shift)
  
  if (is.null(formula_base)) {
    formula_base <- as.formula(paste(dep_var, "~", paste(predictors, collapse = " + ")))
  }
  
  num_candidates <- GA_control$popSize
  candidates <- vector("list", num_candidates)
  set.seed(123)
  for (i in 1:num_candidates) {
    candidate <- list(
      model_approach = sample(genetic_map$model_approach, 1),
      family = sample(genetic_map$family, 1),
      outlier_removal = sample(genetic_map$outlier_removal, 1),
      ranges_of_outlier_removal = sample(genetic_map$ranges_of_outlier_removal, 1),
      outlier_formula = sample(genetic_map$outlier_formula, 1),
      zi_formula = sample(genetic_map$zi_formula, 1),
      normalize_approach = sample(genetic_map$normalize_approach, 1),
      shift_addition = sample(genetic_map$shift_addition, 1),
      shift_ranges = generate_shift_ranges(step_for_shift),
      link = sample(genetic_map$link, 1),
      outlier_removal_approach = sample(genetic_map$outlier_removal_approach, 1),
      zero_inflation_approach = sample(genetic_map$zero_inflation_approach, 1),
      regularization = sample(genetic_map$regularization, 1),
      feature_selection = sample(genetic_map$feature_selection, 1),
      GA_stopping_criterion = sample(genetic_map$GA_stopping_criterion, 1),
      population_initialization = sample(genetic_map$population_initialization, 1),
      cv_folds = sample(genetic_map$cv_folds, 1)
    )
    candidates[[i]] <- candidate
  }
  
  # Use pbapply for candidate evaluation with a progress bar
  results <- pbapply::pblapply(seq_along(candidates), function(i) {
    candidate <- candidates[[i]]
    # Adjust dependent variable for candidate's family requirements
    adj_data <- adjust_dep_var_for_family(data, dep_var, candidate$family)
    
    model <- tryCatch({
      build_model(candidate, adj_data, formula_base, random_effects)
    }, error = function(e) {
      message("Error in build_model (candidate ", i, "): ", e$message)
      NULL
    })
    
    if (!is.null(model)) {
      eval_result <- tryCatch({
        evaluate_fitness(model, relaxation_term = relaxation_term)
      }, error = function(e) {
        message("Error in evaluate_fitness (candidate ", i, "): ", e$message)
        list(fitness = -Inf, p_values = NA, prod_p = NA, sim = NULL)
      })
      list(candidate = candidate, 
           fitness = eval_result$fitness, 
           p_values = eval_result$p_values, 
           prod_p = eval_result$prod_p,
           model = model,
           sim = eval_result$sim,
           outlier_threshold = candidate$ranges_of_outlier_removal)
    } else {
      list(candidate = candidate, fitness = -Inf, p_values = NA, prod_p = NA,
           model = NULL, sim = NULL, outlier_threshold = candidate$ranges_of_outlier_removal)
    }
  })
  
  # Build diagnostic table from candidates that pass the fitness threshold
  diagnostic_list <- lapply(seq_along(results), function(i) {
    res_item <- results[[i]]
    if (!is.null(res_item$model) && !is.na(res_item$fitness) && res_item$fitness != -Inf && res_item$prod_p != 0) {
      data.frame(
        Candidate = paste0("cand_", i),
        ModelApproach = res_item$candidate$model_approach,
        Family = res_item$candidate$family,
        OutlierRemoval = res_item$candidate$outlier_removal,
        RangesOutlierRemoval = res_item$candidate$ranges_of_outlier_removal,
        OutlierFormula = res_item$candidate$outlier_formula,
        ZI_Formula = res_item$candidate$zi_formula,
        NormalizeApproach = res_item$candidate$normalize_approach,
        ShiftAddition = res_item$candidate$shift_addition,
        Link = res_item$candidate$link,
        OutlierRemovalApproach = res_item$candidate$outlier_removal_approach,
        ZeroInflationApproach = res_item$candidate$zero_inflation_approach,
        Regularization = res_item$candidate$regularization,
        FeatureSelection = res_item$candidate$feature_selection,
        GA_StoppingCriterion = res_item$candidate$GA_stopping_criterion,
        PopulationInitialization = res_item$candidate$population_initialization,
        CV_Folds = res_item$candidate$cv_folds,
        Fitness = res_item$fitness,
        ProdP = res_item$prod_p,
        KS = res_item$p_values["ks"],
        Dispersion = res_item$p_values["dispersion"],
        Outlier = res_item$p_values["outlier"],
        Hetero = res_item$p_values["hetero"],
        OutlierThreshold = res_item$outlier_threshold,
        stringsAsFactors = FALSE
      )
    } else {
      NULL
    }
  })
  
  diagnostic_table <- do.call(rbind, diagnostic_list)
  if (!is.null(diagnostic_table)) rownames(diagnostic_table) <- diagnostic_table$Candidate
  
  valid_indices <- which(!is.na(diagnostic_table$Fitness) & diagnostic_table$Fitness != -Inf &
                           diagnostic_table$ProdP != 0)
  
  if (length(valid_indices) > 0) {
    valid_results <- results[valid_indices]
    valid_fitness <- diagnostic_table$Fitness[valid_indices]
    top_order <- order(valid_fitness, decreasing = TRUE)
    top_models <- valid_results[top_order]
    best_model_candidate <- valid_results[[which.max(valid_fitness)]]
    best_model <- best_model_candidate$model
    best_model_sim <- best_model_candidate$sim
    passed_models <- valid_results
  } else {
    top_models <- list()
    best_model <- NULL
    best_model_sim <- NULL
    passed_models <- list()
  }
  
  return(list(all_results = results,
              diagnostic_table = diagnostic_table,
              best_model = best_model,
              best_model_sim = best_model_sim,
              top_models = top_models,
              passed_models = passed_models,
              model_list = results,
              candidates = candidates))
}

#' Demonstration of GA-Based Model Selection
#'
#' Demonstrates the full workflow using the built-in mtcars dataset.
#' Creates a grouping factor for random effects, uses a restricted genetic map,
#' and prints the diagnostic table, best model details, its DHARMa simulation object,
#' and the full list of candidate evaluations.
#'
#' @return Invisibly returns the result list from run_model_selection.
#' @export
demo_model_selection <- function() {
  data(mtcars)
  mtcars$group <- factor(mtcars$cyl)  # Create a grouping factor for random effects
  
  # Define a restricted genetic map for demonstration purposes
  genetic_map <- list(
    model_approach = c("lmer", "glmer", "glmer.nb", "glmmTMB"),
    family = c("gaussian", "binomial", "poisson", "Gamma"),
    outlier_removal = c(TRUE, FALSE),
    ranges_of_outlier_removal = c("<1%", "1%-2.5%", "2.5%-5%", "5%-10%", "10%-15%"),
    outlier_formula = c("~1", "parameter_by_parameter", "all_interactions"),
    zi_formula = c("~1", "single_parameter", "paired_interactions", "full_interactions"),
    normalize_approach = c("no_transform", "arcsinh_x", "boxcox", "log_x", 
                           "double_reverse_log", "sqrt_x", "yeojohnson", "Lambert"),
    shift_addition = c(TRUE, FALSE),
    link = c("log", "logit", "probit", "inverse", "cloglog", "identity", "sqrt"),
    outlier_removal_approach = c("outlier_first", "normalization_first"),
    zero_inflation_approach = c("lmer_nonzero_glmer_binomial", "glmmTMB_zi_formula"),
    regularization = c("none", "L1", "L2", "elastic_net"),
    feature_selection = c("all_predictors", "stepwise", "penalized"),
    GA_stopping_criterion = c("fixed_iterations", "convergence", "max_runtime"),
    population_initialization = c("random", "preselected"),
    cv_folds = c("5_fold", "10_fold", "LOOCV")
  )
  
  res <- run_model_selection(
    data = mtcars,
    dep_var = "mpg",
    predictors = c("wt", "hp"),
    random_effects = "group",
    formula_base = as.formula("mpg ~ wt + hp"),
    genetic_map = genetic_map,
    GA_control = list(popSize = 20, maxiter = 50, run = 25, parallel = TRUE),
    relaxation_term = 0,
    step_for_shift = 0.01,
    missing_values_method = "discard",
    normalize_method = "boxcox"
  )
  
  message("Diagnostic Table:")
  print(res$diagnostic_table)
  
  if (!is.null(res$best_model)) {
    message("\nBest Model:")
    print(res$best_model)
    message("\nDHARMa Simulation Object for Best Model:")
    print(res$best_model_sim)
  } else {
    message("\nNo valid best model was found.")
  }
  
  message("\nFull Candidate Evaluation List (model_list):")
  str(res$model_list)
  
  invisible(res)
}
